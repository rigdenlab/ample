'''Utility file containing functions for ensembling

This structure of the ensembling modules is dictated by the need to be able to pickle
and unpickle the results dictionary. As such all objects need to have a qualified path
(e.g. ensembler.Ensemble ) - otherwise, the module is taken as main, so that when
the results are unpickled, it will look for Ensemble as an attribute of the main ample
script.
'''

import collections
import glob
import iotbx.pdb
import logging
import os
import shutil
import sys

from ample.ensembler import abinitio
from ample.ensembler import homologs
from ample.ensembler import single_model
from ample.util import exit_util
from ample.util import pdb_edit
from ample.util import printTable

__author__ = "Jens Thomas and Felix Simkovic"
__date__ = "18.04.2013"

LOGGER = logging.getLogger(__name__)

def ensembler_factory(amoptd):
    """Return an ensembler object for the required ensembles
    
    :returns: instantiated ensembler object
    """
    if amoptd['homologs']: 
        ensemble_module = homologs
    elif amoptd['single_model_mode']:
        ensemble_module = single_model
    else: 
        ensemble_module = abinitio

    # Paths to ensembling directories
    ensembles_directory = os.path.join(amoptd['work_dir'], 'ensembles')
    amoptd['ensembles_directory'] = ensembles_directory
    work_dir = os.path.join(amoptd['work_dir'], 'ensemble_workdir')

    return ensemble_module.Ensembler(
                                    ensembles_directory=ensembles_directory,
                                    ensemble_max_models=amoptd['ensemble_max_models'],
                                    nproc=amoptd['nproc'],
                                    work_dir=work_dir,
                                    # Executables
                                    fast_protein_cluster_exe=amoptd['fast_protein_cluster_exe'],
                                    gesamt_exe=amoptd['gesamt_exe'],
                                    lsqkab_exe=amoptd['lsqkab_exe'],
                                    maxcluster_exe=amoptd['maxcluster_exe'],
                                    scwrl_exe=amoptd['scwrl_exe'],
                                    spicker_exe=amoptd['spicker_exe'],
                                    theseus_exe=amoptd['theseus_exe'],
                                    )

def cluster_script(amoptd, python_path="ccp4-python"):
    """Create the script for ensembling on a cluster"""
    # write out script
    work_dir = amoptd['work_dir']
    script_path = os.path.join(work_dir, "submit_ensemble.sh")
    with open(script_path, "w") as job_script:
        job_script.write("#!/bin/sh\n")
        job_script.write("ccp4-python -m ample.ensembler -restart_pkl {0}\n".format(amoptd['results_path']))

    # Make executable
    os.chmod(script_path, 0o777)
    return script_path

def create_ensembles(amoptd):
    """Create the ensembles using the values in the amoptd dictionary"""
    
    # Create instance of the ensembler
    ensembler = ensembler_factory(amoptd)

    ############################################################################
    # For a single model we don't need to use glob 
    if not (amoptd['single_model'] or amoptd['models_dir']):
        msg = 'AMPLE ensembler needs either a single_model or a models_dir argument'
        exit_util.exit_error(msg, sys.exc_info()[2])
        if amoptd['single_model'] and not os.path.isfile(amoptd['single_model']):
            msg = 'Cannot find single_model pdb: {0}'.format(amoptd['single_model'])
            exit_util.exit_error(msg, sys.exc_info()[2])
        elif amoptd['models_dir'] and not os.path.isfile(amoptd['models_dir']):
            msg = 'Cannot find models_dir: {0}'.format(amoptd['models_dir'])
            exit_util.exit_error(msg, sys.exc_info()[2])
            
    models = list([amoptd['single_model']]) if amoptd['single_model_mode'] else \
        glob.glob(os.path.join(amoptd['models_dir'], "*.pdb"))

    #if amoptd['cluster_method'] == 'spicker_tmscore':
    #    models = reorder_models(models, amoptd['score_matrix_file_list'])
    #    if not (os.path.isfile(amoptd['score_matrix']) and os.path.isfile(amoptd['score_matrix_file_list'])):
    #        raise RuntimeError("spicker_tmscore needs a score_matrix and score_matrix_file_list")
    #    ensembler.score_matrix = amoptd['score_matrix']

    # Run ensemble creation
    ensembles = ensembler.generate_ensembles_from_amoptd(models, amoptd)

    ############################################################################
    amoptd['ensembles'] = ensembles
    amoptd['ensembles_data'] = ensembler.ensembles_data
    amoptd['truncation_levels'] = ensembler.truncation_levels
    amoptd['truncation_variances'] = ensembler.truncation_variances
    amoptd['truncation_nresidues'] = ensembler.truncation_nresidues

    # We need to let the main process know that we have succeeded as this module could be run on a cluster node with no link
    # to the parent process, so we create a file here indicating that we got this far and didn't die from an exception
    with open(amoptd['ensemble_ok'],'w') as f: f.write('ok\n')

    # Delete all intermediate files if we're purging
    if amoptd['purge']: shutil.rmtree(ensembler.work_dir)
    return

################################################################################
#
# Functions below here are primarily to manipulate, summarise or handle
# existing ensembling data. 
#
################################################################################

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
        exit_util.exit_error(msg)

    LOGGER.info("Importing ensembles from directory: {0}".format(amoptd['ensembles']))

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

def set_phaser_rms_from_subcluster_score(optd):
    """Set the phaser_rms for each ensemble based on its score."""
    if optd['ensemble_options'] is None: optd['ensemble_options'] = {}
    for ensemble_data in optd['ensembles_data']:
        name = ensemble_data['name']
        if name not in optd['ensemble_options']: optd['ensemble_options'][name] = {}
        optd['ensemble_options'][name]['phaser_rms'] = ensemble_data['subcluster_score']
        
    return

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
    tmp_data = collections.defaultdict(list)
    if "cluster_num" in keys_to_sort:
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['cluster_num']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys(), key=int)

    elif "truncation_score_key" in keys_to_sort:
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['truncation_score_key']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())

    else:
        for ensemble in ensembles_zipped:
            tmp_data["all"].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())

    ############################################################################
    #
    # Iterate through clusters and group based on truncation level
    #
    ############################################################################
    ensembles_zipped_ordered = []

    for sbin in iterator_keys:
        low, mid, high = [], [], []
        for ensemble in _sort_ensembles_parameters(tmp_data[sbin], keys_to_sort):
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
