'''
Created on Apr 18, 2013

This structure of the ensembling modules is dictated by the need to be able to pickle
and unpickle the results dictionary. As such all objects need to have a qualified path
(e.g. ample_ensemble.Ensemble ) - otherwise, the module is taken as main, so that when
the results are unpickled, it will look for Ensemble as an attribute of the main ample
script.

@author: jmht
'''

import cPickle
import glob
import itertools
import logging
import os
import random
import shutil
import sys
import unittest

# our imports
import ample_ensemble
import ample_exit
import ample_util
import iotbx.pdb
import pdb_edit
import printTable

def cluster_script(amoptd, python_path="ccp4-python"):
    """Create the script for ensembling on a cluster"""
    # write out script
    work_dir = amoptd['work_dir']
    script_path = os.path.join(work_dir, "submit_ensemble.sh")
    with open(script_path, "w") as job_script:
        job_script.write("#!/bin/sh\n")
        # Find path to this directory to get path to python ensemble.py script
        pydir = os.path.abspath(os.path.dirname(__file__))
        ensemble_script = os.path.join(pydir, "ensemble.py")
        job_script.write("{0} {1} {2} {3}\n".format(python_path, "-u", ensemble_script, amoptd['results_path']))

    # Make executable
    os.chmod(script_path, 0o777)
    return script_path

def create_ensembles(amoptd):
    """Create the ensembles using the values in the amoptd dictionary"""

    ensembler = ample_ensemble.Ensembler()
    
    ensembler.theseus_exe = amoptd['theseus_exe']
    if amoptd['subcluster_program']:
        ensembler.subcluster_program = amoptd['subcluster_program']
        if amoptd['subcluster_program'] == 'maxcluster':
            ensembler.subcluster_exe = amoptd['maxcluster_exe']
        elif amoptd['subcluster_program'] == 'lsqkab':
            ensembler.subcluster_exe = ensembler.lsqkab_exe
        else:
            raise RuntimeError('Unrecognised subcluster_program: {0}'.format(amoptd['subcluster_program']))
        
    ensembler.max_ensemble_models = amoptd['max_ensemble_models']
    
    if amoptd['cluster_method'] == 'spicker' or \
       amoptd['cluster_method'] == 'spicker_qscore' or \
       amoptd['cluster_method'] == 'spicker_tmscore':
        cluster_exe = amoptd['spicker_exe']
        if amoptd['cluster_method'] == 'spicker_tmscore':
            if not (os.path.isfile(amoptd['score_matrix']) and os.path.isfile(amoptd['score_matrix_file_list'])):
                    raise RuntimeError("spicker_tmscore needs a score_matrix and score_matrix_file_list")
            ensembler.score_matrix = amoptd['score_matrix']
    elif amoptd['cluster_method'] == 'fast_protein_cluster':
        cluster_exe = amoptd['fast_protein_cluster_exe']
    elif amoptd['cluster_method'] == 'import' or \
         amoptd['cluster_method'] == 'random':
        cluster_exe = None
    else:
        raise RuntimeError, "create_ensembles - unrecognised cluster_method: {0}".format(amoptd['cluster_method'])
    ensembler.scwrl_exe = amoptd['scwrl_exe']
    ensembler.gesamt_exe = amoptd['gesamt_exe']
        
    ensembles_directory = os.path.join(amoptd['work_dir'], 'ensembles')
    if not os.path.isdir(ensembles_directory): os.mkdir(ensembles_directory)
    amoptd['ensembles_directory'] = ensembles_directory
    
    work_dir = os.path.join(amoptd['work_dir'], 'ensemble_workdir')
    os.mkdir(work_dir)
    os.chdir(work_dir)

    models = glob.glob(os.path.join(amoptd['models_dir'], "*.pdb"))
    if amoptd['cluster_method'] == 'spicker_tmscore':   
        models = reorder_models(models, amoptd['score_matrix_file_list'])
    
    if amoptd['homologs']:
        ensembles = ensembler.generate_ensembles_homologs(models,
                                                          percent_truncation=amoptd['percent'],
                                                          truncation_method=amoptd['truncation_method'],
                                                          ensembles_directory=ensembles_directory,
                                                          side_chain_treatments=amoptd['side_chain_treatments'],
                                                          alignment_file=amoptd['alignment_file'],
                                                          work_dir=work_dir,
                                                          nproc=amoptd['nproc'],
                                                          homolog_aligner=amoptd['homolog_aligner'],
                                                          mustang_exe=amoptd['mustang_exe'],
                                                          gesamt_exe=amoptd['gesamt_exe'],
                                                          use_scwrl=amoptd['use_scwrl'])
    else:
        ensembles = ensembler.generate_ensembles(models,
                                                 cluster_method=amoptd['cluster_method'],
                                                 cluster_exe=cluster_exe,
                                                 num_clusters=amoptd['num_clusters'],
                                                 cluster_dir=amoptd['cluster_dir'],
                                                 percent_truncation=amoptd['percent'],
                                                 truncation_method=amoptd['truncation_method'],
                                                 truncation_pruning=amoptd['truncation_pruning'],
                                                 ensembles_directory=ensembles_directory,
                                                 side_chain_treatments=amoptd['side_chain_treatments'],
                                                 use_scwrl=amoptd['use_scwrl'],
                                                 work_dir=work_dir,
                                                 nproc=amoptd['nproc'])
    
    amoptd['ensembles'] = ensembles
    amoptd['ensembles_data'] = ensembler.ensembles_data
    amoptd['truncation_levels'] = ensembler.truncation_levels
    amoptd['truncation_variances'] = ensembler.truncation_variances
    amoptd['truncation_nresidues'] = ensembler.truncation_nresidues
    
    # We need to let the main process know that we have succeeded as this module could be run on a cluster node with no link
    # to the parent process, so we create a file here indicating that we got this far and didn't die from an exception
    with open(amoptd['ensemble_ok'],'w') as f: f.write('ok\n')
    
    # Delete all intermediate files if we're purging
    if amoptd['purge']: shutil.rmtree(work_dir)
    return

def collate_cluster_data(ensembles_data):
    clusters = {}  # Loop through all ensemble data objects and build up a data tree
    cluster_method = None
    truncation_method = None
    percent_truncation = None
    side_chain_treatments = []
    for e in ensembles_data:
        if not cluster_method:
            cluster_method = e['cluster_method']
            percent_truncation = e['percent_truncation']
            truncation_method = e['truncation_method']
            # num_clusters = e['num_clusters']
        cnum = e['cluster_num']
        if cnum not in clusters:
            clusters[cnum] = {}
            clusters[cnum]['cluster_centroid'] = e['cluster_centroid']
            clusters[cnum]['cluster_num_models'] = e['cluster_num_models']
            clusters[cnum]['tlevels'] = {}
        tlvl = e['truncation_level']
        if tlvl not in clusters[cnum]['tlevels']:
            clusters[cnum]['tlevels'][tlvl] = {}
            clusters[cnum]['tlevels'][tlvl]['truncation_variance'] = e['truncation_variance']
            clusters[cnum]['tlevels'][tlvl]['num_residues'] = e['num_residues']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'] = {}
        srt = e['subcluster_radius_threshold']
        if srt not in clusters[cnum]['tlevels'][tlvl]['radius_thresholds']:
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt] = {}
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['num_models'] = e['subcluster_num_models']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'] = {}
        sct = e['side_chain_treatment']
        if sct not in side_chain_treatments: side_chain_treatments.append(sct)
        if sct not in clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct']:
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct] = {}
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['name'] = e['name']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['num_atoms'] = e['ensemble_num_atoms']
    
    return clusters, cluster_method, truncation_method, percent_truncation, side_chain_treatments

def cluster_table_data(clusters, cluster_num, side_chain_treatments):
    # FIX TO IGNORE side_chain_treatments
    # tdata = [("Name", "Truncation Level", u"Variance Threshold (\u212B^2)", "No. Residues", u"Radius Threshold (\u212B)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
    tdata = [("Name", "Cluster", "Truncation Level", "Variance Threshold (A^2)", "No. Residues", "Radius Threshold (A)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
    for tl in sorted(clusters[cluster_num]['tlevels']):
        tvar = clusters[cluster_num]['tlevels'][tl]['truncation_variance']
        nresidues = clusters[cluster_num]['tlevels'][tl]['num_residues']
        for i, rt in enumerate(sorted(clusters[cluster_num]['tlevels'][tl]['radius_thresholds'])):
            nmodels = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['num_models']
            side_chain_treatments = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'].keys()
            for j, sct in enumerate(side_chain_treatments):
                name = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['name']
                num_atoms = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['num_atoms']
                if i == 0 and j == 0:  # change of radius
                    tdata.append((name, cluster_num, tl, tvar, nresidues, rt, nmodels, num_atoms, sct))
                elif i > 0 and j == 0:  # change of side_chain
                    tdata.append((name, "", "", "", "", rt, nmodels, num_atoms, sct))
                else:
                    tdata.append((name, "", "", "", "", "", "", num_atoms, sct))
    return tdata

def ensemble_summary(ensembles_data):
    """Print a summary of the ensembling process"""

    clusters, cluster_method, truncation_method, percent_truncation, side_chain_treatments = collate_cluster_data(ensembles_data)
    num_clusters = len(clusters)
    
    tableFormat = printTable.Table()
    rstr = "\n"
    rstr += "Ensemble Results\n"
    rstr += "----------------\n\n"
    rstr += "Cluster method: {0}\n".format(cluster_method)
    rstr += "Truncation method: {0}\n".format(truncation_method)
    rstr += "Percent truncation: {0}\n".format(percent_truncation)
    rstr += "Number of clusters: {0}\n".format(num_clusters)
    
    for cluster_num in sorted(clusters.keys()):
        rstr += "\n"
        rstr += "Cluster {0}\n".format(cluster_num)
        rstr += "Number of models: {0}\n".format(clusters[cluster_num]['cluster_num_models'])
        rstr += "Cluster centroid: {0}\n".format(clusters[cluster_num]['cluster_centroid'])
        rstr += "\n"
        tdata = cluster_table_data(clusters, cluster_num, side_chain_treatments)
        rstr += tableFormat.pprint_table(tdata)        
    
    rstr += "\nGenerated {0} ensembles\n\n".format(len(ensembles_data))
    
    return rstr

def import_ensembles(amoptd):
    if not pdb_edit.check_pdb_directory(amoptd['ensembles'], single=False):
        msg = "Cannot import ensembles from the directory: {0}".format(amoptd['ensembles'])
        ample_exit.exit_error(msg)

    msg = "Importing ensembles from directory: {0}".format(amoptd['ensembles'])
    logging.info(msg)
    ensembles = glob.glob(os.path.join(amoptd['ensembles'], '*.pdb'))
    amoptd['ensembles'] = ensembles
    
    # get the data on the ensemble
    ensembles_data = []
    for e in ensembles:
        d = {}
        d['name'] = os.path.splitext(os.path.basename(e))[0]
        d['ensemble_pdb'] = e
        
        # Get data on the models
        hierarchy = iotbx.pdb.pdb_input(file_name=e).construct_hierarchy()
        d['subcluster_num_models'] = len(hierarchy.models())
        d['num_residues'] = len(hierarchy.models()[0].chains()[0].residue_groups())
        d['ensemble_num_atoms'] = len(hierarchy.models()[0].atoms())
        
        ensembles_data.append(d)
        
    amoptd['ensembles_data'] = ensembles_data
        
    return ensembles

def reorder_models(models, ordered_list_file):
    """Reorder the list of models from a list of models (possibly in a different directory"""
    with open(ordered_list_file) as f:
        ordered_list = [ line.strip() for line in f if line.strip() ]
    mdir = os.path.dirname(models[0])
    model_names = set([ os.path.basename(m) for m in models ])
    reordered_models = []
    for tm_model in ordered_list:
        name = os.path.basename(tm_model)
        assert name in model_names,"TM model {0} is not in models!".format(name)
        reordered_models.append(os.path.join(mdir, name))
    return reordered_models

def sort_ensembles(ensemble_pdbs, ensembles_data=None, prioritise=True):
    """Sort AMPLE ensemble data.
    
    :ensemble_pdbs: list of ensembles
    :ensembles_data: data found in amoptd['ensembles_data'] entry
    :prioritise: default True; favour ensembles in truncation level 20-50%
    
    :returns: list of sorted ensembles pdbs
    """
    if len(ensemble_pdbs) < 1: raise RuntimeError("Not enough ensembles")

    # Sort with out method only if we have data which contains our generated data
    if ensembles_data:
        ensemble_pdbs_sorted = _sort_ensembles(ensemble_pdbs, ensembles_data,
                                               prioritise)    
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)
    
    return ensemble_pdbs_sorted

def _sort_ensembles(ensemble_pdbs, ensembles_data, prioritise):
    """Sort ensembles based on
        1) Cluster
        2) Truncation level - favoured region of 20-50% first
        3) Subcluster radius threshold
        4) Side chain treatment
    """
    assert len(ensemble_pdbs) == len(ensembles_data), "Unequal ensembles data for sorting"
    
    # Determine which keys we can sort our data with
    def _get_keys(ensemble):
        _criteria = ['cluster_num', 
                     'truncation_score_key', 
                     'truncation_level', 
                     'subcluster_radius_threshold',
                     'side_chain_treatment']
        
        for crit in _criteria:
            if ensemble.has_key(crit) and ensemble[crit]:
                yield crit
    
    # Keys we want to sort data by - differs for different source of ensembles
    keys_to_sort = list(_get_keys(ensembles_data[0])) 
    
    if keys_to_sort:
        # Zip the data so order remains identical between pdbs and data
        ensembles_zipped = zip(ensemble_pdbs, ensembles_data)
    
        if "truncation_level" in keys_to_sort and prioritise:
            ensembles_zipped_ordered = _sort_ensembles_prioritise(ensembles_zipped, keys_to_sort)
        else:
            ensembles_zipped_ordered = _sort_ensembles_parameters(ensembles_zipped, keys_to_sort)
    
        # We need to `unzip` the data to get a list of ordered pdbs
        ensemble_pdbs_sorted, _ = zip(*ensembles_zipped_ordered)
    
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)
    
    return ensemble_pdbs_sorted
        
def _sort_ensembles_prioritise(ensembles_zipped, keys_to_sort):
    """Sort our ensembles based on the data in order of the keys provided"""
    
    ############################################################################
    #
    # Group the ensembles either by cluster ... or by truncation score ...
    # ... otherwise throw them all in one pot
    #
    ############################################################################
    if "cluster_num" in keys_to_sort:
        tmp_data = {ens[1]['cluster_num']: [] for ens in ensembles_zipped}
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['cluster_num']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys(), key=int)
        
    elif "truncation_score_key" in keys_to_sort:
        tmp_data = {ens[1]['truncation_score_key']: [] for ens in ensembles_zipped}
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['truncation_score_key']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())
        
    else:
        tmp_data = {"all": []}
        for ensemble in ensembles_zipped:
            tmp_data["all"].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())
    
    ############################################################################
    #
    # Iterate through clusters and group based on truncation level
    #
    ############################################################################
    ensembles_zipped_ordered = []
    
    for bin in iterator_keys:
        low, mid, high = [], [], []
        for ensemble in _sort_ensembles_parameters(tmp_data[bin], keys_to_sort):
            if ensemble[1]['truncation_level'] > 50: 
                high.append(ensemble)
            elif ensemble[1]['truncation_level'] < 20: 
                low.append(ensemble)
            else: 
                mid.append(ensemble)
        ensembles_zipped_ordered += mid + low + high
    
    return ensembles_zipped_ordered
 
def _sort_ensembles_parameters(ensembles_zipped, keys_to_sort):
    """Tiny wrapper function to sort ensembles data based on keys provided
    """
    def _extract(ens):
        return [ens[1][crit] for crit in keys_to_sort]
    return sorted(ensembles_zipped, key=_extract)


def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testSummary'))
    return suite

class Test(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        cls.theseus_exe = ample_util.find_exe('theseus')
        cls.spicker_exe = ample_util.find_exe('spicker')
        cls.maxcluster_exe = ample_util.find_exe('maxcluster')

        return
    
    def testSummary(self):

        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = ample_ensemble.Ensembler()

        work_dir = os.path.join(os.getcwd(), "summary")
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")

        num_clusters = 1
        cluster_method = 'spicker'
        percent_truncation = 20
        truncation_method = "percent"
        ensembler.generate_ensembles(models,
                                     cluster_method=cluster_method,
                                     cluster_exe=self.spicker_exe,
                                     num_clusters=num_clusters,
                                     percent_truncation=percent_truncation,
                                     truncation_method=truncation_method,
                                     work_dir=work_dir)
        
        print ensemble_summary(ensembler.ensembles_data)
        shutil.rmtree(work_dir)
        
        return

    def test_sortEnsembles1(self):
        """Test sorting for ab initio data"""
        
        clusters = [1, 2, 3]
        tlevels = [19, 20, 50, 80, 100]
        radii = [1, 2, 3]
        sidechains = ["allatom", "polyAla", "reliable"]
        
        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(clusters)
        random.shuffle(tlevels)
        random.shuffle(radii)
        random.shuffle(sidechains)
        
        ensemble_data = []
        ensemble_pdbs = []
        
        # Generate ensembles to sort 
        for comb in list(itertools.product(clusters, tlevels, radii, sidechains)):
            
            name = "c{cluster}_tl{truncation}_r{radius}_{schain}".format(cluster=comb[0], 
                                                                         truncation=comb[1], 
                                                                         radius=comb[2], 
                                                                         schain=comb[3])
            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'cluster_num': comb[0],
                        'truncation_level': comb[1],
                        'subcluster_radius_threshold': comb[2],
                        'side_chain_treatment': comb[3],
            }
            
            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)
        
        ensemble_pdb_sorted_1 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        ensemble_pdb_sorted_2 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        ensemble_pdb_sorted_3 = sort_ensembles(ensemble_pdbs)
        
        self.assertEqual("/foo/bar/c1_tl20_r1_allatom.pdb", ensemble_pdb_sorted_1[0])
        self.assertEqual("/foo/bar/c1_tl19_r1_allatom.pdb", ensemble_pdb_sorted_1[18])
        self.assertEqual("/foo/bar/c1_tl100_r3_reliable.pdb", ensemble_pdb_sorted_1[44])
        self.assertEqual("/foo/bar/c2_tl20_r1_allatom.pdb", ensemble_pdb_sorted_1[45])
        self.assertEqual("/foo/bar/c2_tl19_r1_allatom.pdb", ensemble_pdb_sorted_1[63])
        self.assertEqual("/foo/bar/c2_tl100_r3_reliable.pdb", ensemble_pdb_sorted_1[89])
        self.assertEqual("/foo/bar/c3_tl20_r1_allatom.pdb", ensemble_pdb_sorted_1[90])    
        self.assertEqual("/foo/bar/c3_tl19_r1_allatom.pdb", ensemble_pdb_sorted_1[108])
        self.assertEqual("/foo/bar/c3_tl100_r3_reliable.pdb", ensemble_pdb_sorted_1[-1])
        
        self.assertEqual("/foo/bar/c1_tl19_r1_allatom.pdb", ensemble_pdb_sorted_2[0])
        self.assertEqual("/foo/bar/c3_tl100_r3_reliable.pdb", ensemble_pdb_sorted_2[-1])    
        self.assertEqual("/foo/bar/c2_tl50_r2_polyAla.pdb", ensemble_pdb_sorted_2[67])
        
        self.assertEqual("/foo/bar/c1_tl100_r1_allatom.pdb", ensemble_pdb_sorted_3[0])
        self.assertEqual("/foo/bar/c3_tl80_r3_reliable.pdb", ensemble_pdb_sorted_3[-1])
        
        return
    
    def test_sortEnsembles2(self):
        """Test sorting for homolog data"""
        
        tlevels = [19, 20, 50, 80, 100]
        sidechains = ["allatom", "polyAla", "reliable"]
        
        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(tlevels)
        random.shuffle(sidechains)
        
        ensemble_data = []
        ensemble_pdbs = []
        
        # Generate ensembles to sort 
        for comb in list(itertools.product(tlevels, sidechains)):
            
            name = "e{truncation}_{schain}".format(truncation=comb[0], 
                                                   schain=comb[1])
            
            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_level': comb[0],
                        'side_chain_treatment': comb[1],
            }
            
            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)
        
        ensemble_pdb_sorted_1 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        ensemble_pdb_sorted_2 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        ensemble_pdb_sorted_3 = sort_ensembles(ensemble_pdbs)
        
        self.assertEqual("/foo/bar/e20_allatom.pdb", ensemble_pdb_sorted_1[0])
        self.assertEqual("/foo/bar/e19_polyAla.pdb", ensemble_pdb_sorted_1[len(ensemble_pdb_sorted_1)/2])
        self.assertEqual("/foo/bar/e100_reliable.pdb", ensemble_pdb_sorted_1[-1])
        
        self.assertEqual("/foo/bar/e19_allatom.pdb", ensemble_pdb_sorted_2[0])
        self.assertEqual("/foo/bar/e100_reliable.pdb", ensemble_pdb_sorted_2[-1])    
        self.assertEqual("/foo/bar/e50_polyAla.pdb", ensemble_pdb_sorted_2[len(ensemble_pdb_sorted_1)/2])
        
        self.assertEqual("/foo/bar/e100_allatom.pdb", ensemble_pdb_sorted_3[0])
        self.assertEqual("/foo/bar/e80_reliable.pdb", ensemble_pdb_sorted_3[-1])
        
        return
    
    def test_sortEnsembles3(self):
        '''Test sorting for single model data'''
        
        score_keys = ["ntv", "bbc", "kicker"]
        tlevels = [19, 20, 50, 80, 100]
        sidechains = ["allatom", "polyAla", "reliable"]
        
        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(score_keys)
        random.shuffle(tlevels)
        random.shuffle(sidechains)
        
        ensemble_data = []
        ensemble_pdbs = []
        
        # Generate ensembles to sort 
        for comb in list(itertools.product(score_keys, tlevels, sidechains)):
            
            name = "{score_key}_tl{truncation}_{schain}".format(score_key=comb[0], 
                                                                truncation=comb[1], 
                                                                schain=comb[2])
            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_score_key': comb[0],
                        'truncation_level': comb[1],
                        'side_chain_treatment': comb[2],
            }
            
            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)
        
        ensemble_pdb_sorted_1 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        ensemble_pdb_sorted_2 = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        ensemble_pdb_sorted_3 = sort_ensembles(ensemble_pdbs)
        
        self.assertEqual("/foo/bar/bbc_tl20_allatom.pdb", ensemble_pdb_sorted_1[0])
        self.assertEqual("/foo/bar/bbc_tl19_allatom.pdb", ensemble_pdb_sorted_1[6])
        self.assertEqual("/foo/bar/bbc_tl100_reliable.pdb", ensemble_pdb_sorted_1[14])
        self.assertEqual("/foo/bar/kicker_tl20_allatom.pdb", ensemble_pdb_sorted_1[15])
        self.assertEqual("/foo/bar/kicker_tl19_allatom.pdb", ensemble_pdb_sorted_1[21])
        self.assertEqual("/foo/bar/kicker_tl100_reliable.pdb", ensemble_pdb_sorted_1[29])
        self.assertEqual("/foo/bar/ntv_tl20_allatom.pdb", ensemble_pdb_sorted_1[30])    
        self.assertEqual("/foo/bar/ntv_tl19_allatom.pdb", ensemble_pdb_sorted_1[36])
        self.assertEqual("/foo/bar/ntv_tl100_reliable.pdb", ensemble_pdb_sorted_1[-1])
        
        self.assertEqual("/foo/bar/bbc_tl19_allatom.pdb", ensemble_pdb_sorted_2[0])
        self.assertEqual("/foo/bar/ntv_tl100_reliable.pdb", ensemble_pdb_sorted_2[-1])    
        self.assertEqual("/foo/bar/kicker_tl50_polyAla.pdb", ensemble_pdb_sorted_2[len(ensemble_pdb_sorted_2)/2])
        
        self.assertEqual("/foo/bar/bbc_tl100_allatom.pdb", ensemble_pdb_sorted_3[0])
        self.assertEqual("/foo/bar/ntv_tl80_reliable.pdb", ensemble_pdb_sorted_3[-1])
        
        return
    
if __name__ == "__main__":

    # This runs the ensembling starting from a pickled file containing an amopt dictionary.
    # - used when submitting the modelling jobs to a cluster

    if len(sys.argv) != 2 or not os.path.isfile(sys.argv[1]):
        print "ensemble script requires the path to a pickled amopt dictionary!"
        sys.exit(1)

    # Get the amopt dictionary
    with open(sys.argv[1], "r") as f: amoptd = cPickle.load(f)

    # if os.path.abspath(fpath) != os.path.abspath(amoptd['results_path']):
    #    print "results_path must match the path to the pickle file"
    #    sys.exit(1)

    # Set up logging - could append to an existing log?
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fl = logging.FileHandler(os.path.join(amoptd['work_dir'],"ensemble.log"))
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)

    # Create the ensembles & save them
    create_ensembles(amoptd)
    ample_util.saveAmoptd(amoptd)
