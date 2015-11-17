'''
Created on Apr 18, 2013

@author: jmht
'''

import collections
import copy
import glob
import logging
import os
import random
import re
import shutil
import sys
import unittest

# our imports
import ample_sequence
import ample_util
import fast_protein_cluster
import pdb_edit
import spicker
import subcluster
import theseus

POLYALA = 'polyAla'
RELIABLE = 'reliable'
ALLATOM = 'allatom'
UNMODIFIED = 'unmod'
SIDE_CHAIN_TREATMENTS = [POLYALA, RELIABLE, ALLATOM]

def align_mustang(models, mustang_exe=None, work_dir=None):
    if not ample_util.is_exe(mustang_exe):
        raise RuntimeError, "Cannot find mustang executable: {0}".format(mustang_exe)
    
    if not work_dir: work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    os.chdir(work_dir)

    logfile = os.path.join(work_dir, 'mustang.log')
    basename = 'mustang'
    cmd = [mustang_exe, '-F', 'fasta', '-o', basename, '-i' ] + models
    rtn = ample_util.run_command(cmd, logfile=logfile, directory=work_dir)
    if not rtn == 0:
        raise RuntimeError, "Error running mustang. Check logfile: {0}".format(logfile)
    
    alignment_file = os.path.join(work_dir, basename + ".afasta")
    if not os.path.isfile(alignment_file): raise RuntimeError, "Could not find alignment file: {0} after running mustang!".format(alignment_file)
    return alignment_file

def align_gesamt(models, gesamt_exe=None, work_dir=None):
    if not ample_util.is_exe(gesamt_exe):
        raise RuntimeError, "Cannot find gesamt executable: {0}".format(gesamt_exe)
    
    if not work_dir: work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    os.chdir(work_dir)
    
    # Need to map chain name to pdb
    model2chain = {}
    for m in models:
        seqd = pdb_edit.sequence(m)
        if len(seqd) != 1: raise RuntimeError, "Model {0} does not contain a single chain, got: {1}".format(seqd.keys())
        model2chain[m] = seqd.keys()[0]
    
    basename = 'gesamt'
    logfile = os.path.join(work_dir, 'gesamt.log')
    alignment_file = os.path.join(work_dir, basename + ".afasta")
    
    # Build up command-line
    cmd = [gesamt_exe]
    # We iterate through the models to make sure the order stays the same
    for m in models: cmd += [ m, '-s', model2chain[m] ]
    cmd += ['-o', '{0}.pdb'.format(basename), '-a', alignment_file]
    
    rtn = ample_util.run_command(cmd, logfile=logfile, directory=work_dir)
    if not rtn == 0:
        raise RuntimeError, "Error running gesamt. Check logfile: {0}".format(logfile)
    
    if not os.path.isfile(alignment_file): raise RuntimeError, "Could not find alignment file: {0} after running gesamt!".format(alignment_file)
    return alignment_file
    
def model_core_from_alignment(models, alignment_file, work_dir=None):
    
    if not work_dir: work_dir = os.path.join(os.getcwd(), 'core_models')
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    
    # Read in alignment to get
    align_seq = ample_sequence.Sequence(fasta=alignment_file)
    
    # Check all alignments the same length
    
    # Get pdb names from alignment headers
    seq_names = [ h[1:].strip() for h in align_seq.headers ]
    
    # Need to check if the alignment file is from gesamt, in which case, the names have the
    # chain names in brackets appended
    for i, s in enumerate(seq_names):
        x = re.search("\([a-zA-Z]*\)$", s)
        if x: seq_names[i] = s.replace(x.group(0), "")
    
    # Get array specifying which positions are core. If the positions all align, then there
    # will be a capital letter for the residue. Gaps are signified by "-" and non-structurally-
    # aligned residues by lower-case letters
    GAP = '-'
    # Can't use below as Theseus ignores lower-case letters in the alignment
    # core = [ all([ x in pdb_edit.one2three.keys() for x in t ]) for t in zip(*align_seq.sequences) ]
    core = [ all([ x != GAP for x in t ]) for t in zip(*align_seq.sequences) ]
    
    # For each sequence, get a list of which positions are core
    core_positions = []
    for seq in align_seq.sequences:
        p = []
        count = 0
        for i, pos in enumerate(seq):
            if pos != GAP:
                if core[i]: p.append(count)
                count += 1
        core_positions.append(p)
        
    # Should check lengths of sequences match the length of the aa in the pdbs
        
    # Create dict mapping seq_names to core positions
    core_dict = dict((s, core_positions[i]) for i, s in enumerate(seq_names))
    
    # Cut the models down to core
    core_models = []
    for m in models:
        name = os.path.basename(m)
        pdbout = ample_util.filename_append(m, astr='core', directory=work_dir)
        pdb_edit.select_residues(m, pdbout, tokeep_idx=core_dict[name])
        core_models.append(pdbout)
        
    return core_models

def model_core_from_theseus(models, alignment_file, var_by_res, work_dir=None):
    """Only residues from the first protein are listed in the theseus output, but then not even all of them
    
    We assume the output is based on the original alignment so that where each residue in the first protein lines up with either another residue in one of the other proteins or a gap
    
    SO - we need to go through the theseus data and for each residue that is core find the corresponding residues in the other proteins
    
    We use the resSeq numbers to match the residues across the alignment
    """
    if not work_dir: work_dir = os.path.join(os.getcwd(), 'core_models')
    if not os.path.isdir(work_dir): os.mkdir(work_dir)

    seqalign = ample_sequence.Sequence(fasta=alignment_file)

    # Get the names of the pdb files from the fasta header
    # Format is expected to be: '>1ujb.pdb(A)'
    names = [ h[1:].split('(')[0] for h in seqalign.headers ]
    first = names[0]
    
    # Dictionary mapping model pdb to resSeqs that are core
    model2core = {}
    for n in names: model2core[n] = [] # initialise
    
    # We now need to add the list of resSeqs of the other models to the Sequence object
    for m in models: seqalign.add_pdb_data(m)
    
    # Get list of core resSeqs in the first sequence
    model2core[first] = [ x.resSeq for x in var_by_res if x.core ]
    
    # Now go through the first sequence and get the resSeqs of the corresponding core for the other models
    needle = 0
    for i, resSeq in enumerate(seqalign.resseqs[0]):
        if model2core[first][needle] == resSeq:
            #print i,resSeq, needle, model2core[first][needle]
            # Core residue in first sequence so append the corresponding resSequs for the other proteins
            for j, n in enumerate(names[1:]):
                model2core[n].append(seqalign.resseqs[j+1][i])
            needle += 1
            if needle >= len(model2core[first]): break
    
    work_dir = os.getcwd()
    
    core_models = []
    for m in models:
        name = os.path.basename(m)
        pdbout = ample_util.filename_append(m, astr='core', directory=work_dir)
        print "MKAING ",pdbout
        pdb_edit.select_residues(m, pdbout, tokeep_idx=model2core[name])
        core_models.append(pdbout)
        
    return core_models

def split_sequence(length, percent_interval, min_chunk=3):
    """split a sequence of length into chunks each separated by percent_interval each being at least min_chunk size"""
    
    # How many residues should fit in each bin
    chunk_size = int(round(float(length) * float(percent_interval) / 100.0))
    idxs = [length - 1]
    while True:
        start = idxs[-1] - chunk_size
        remainder = start + 1
        if remainder >= min_chunk:
            idxs.append(start)
        else:
            break
    return idxs

class Ensembler(object):
    """Class to generate ensembles from ab inito models (all models must have same sequence)
    
    """
    def __init__(self):
        
        # For all
        self.work_dir = None  # top directory where everything gets done
        self.theseus_exe = None
        
        # clustering
        self.cluster_method = "spicker"  # the method for initial clustering
        self.cluster_exe = None
        self.num_clusters = 1
        
        # truncation
        self.truncation_method = "percent"
        self.truncation_pruning = "none"
        self.percent_truncation = 5
        self.pruning_strategy = "none"
        # For Felix so we know what truncation levels we used
        self.truncation_levels = None
        self.truncation_variances = None
        self.truncation_nresidues = None
        
        # subclustering
        # self.subcluster_method='FLOATING_RADII'
        self.subcluster_method = 'ORIGINAL'
        self.subcluster_program = "maxcluster"
        self.subcluster_exe = None
        self.subclustering_method = "radius"
        self.subcluster_radius_thresholds = [1, 2, 3]
        self.ensemble_max_models = 30
        
        # ensembles
        self.ensembles_directory = None
        self.ensembles = None
        self.ensembles_data = None
        
        # misc
        self.logger = logging.getLogger()
        
        return
    
    def align_models(self, models, basename=None, work_dir=None, homologs=False):
        run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.align_models(models, basename=basename, homologs=homologs)
        except Exception, e:
            self.logger.critical("Error running theseus: {0}".format(e))
            return False
        return run_theseus.superposed_models
    
    def _calculate_residues_focussed(self, var_by_res):
        """
        The sweet spot for success seems to occur in the interval 5-40 residues.
        Up till now we have always worked in 5% intervals, so 20 truncation levels
        The new strategy is to ensure that always have at least half of the truncations in
        the interval < 40 residues => 10 truncations in 40, so at least 4 residue chunks in this interval.
        
        The strategy is therefore for < 80 residues, just split evenly into 20 chunks.
        
        For > 80 residues, split < 40 into 10 4-residue chunks, and split the interval 40 -> end into
        10 even chunks.
        """
        
        length = len(var_by_res)
        if length <= 80:
            # Just split evenly into 20 chunks
            return self._calculate_residues_percent(var_by_res, 5)
    
        # Get list of residue indices sorted by variance - from least variable to most
        var_by_res.sort(key=lambda x: x.variance, reverse=False)
         
        # Split a 40 - length interval into 10 even chunks.
        llen = 40
        lower_start = split_sequence(llen, 10)
        
        # Split remaining interval into 10 even chunks. We need to add the start sequence as we have
        # removed llen residues
        ulen = length - llen
        upper_start = [ i + llen for i in split_sequence(ulen, 10) ]
        start_indexes = upper_start + lower_start 
        
        # Calculate the percentages for each of these start points
        percentages = [ int(round(float(start + 1) / float(length) * 100)) for start in start_indexes ]
        # print "percentages ", percentages
        truncation_levels = percentages

        # print "var_by_res ",var_by_res
        idxs_all = [ x.idx for x in var_by_res ]
        resseq_all = [ x.resSeq for x in var_by_res ]
        variances = [ x.variance for x in var_by_res ]

        truncation_residue_idxs = [ sorted(idxs_all[:i + 1]) for i in start_indexes ]
        # print "truncation_residue_idxs ",truncation_residue_idxs
        truncation_residues = [ sorted(resseq_all[:i + 1]) for i in start_indexes ]
        # print "truncation_residues ",truncation_residues
        
        # We take the variance of the most variable residue
        truncation_variances = [ variances[i] for i in start_indexes ] 
        # print "truncation_variances ",truncation_variances
        
        return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs

    def _calculate_residues_percent(self, var_by_res, percent_interval):
        """Calculate the list of residues to keep if we are keeping self.percent residues under
        each truncation bin. The threshold is just the threshold of the most variable residue"""
        
        MIN_CHUNK = 3  # We need at least 3 residues for theseus to work
        length = len(var_by_res)
        start_idxs = split_sequence(length, percent_interval, min_chunk=MIN_CHUNK)
        
        # Get list of residue indices sorted by variance - from least to most
        print "GOT ",var_by_res
        print "GOT ",var_by_res[0]
        var_by_res.sort(key=lambda x: x.variance, reverse=False)
         
        # print "var_by_res ",var_by_res
        idxs_all = [ x.idx for x in var_by_res ]
        resseq_all = [ x.resSeq for x in var_by_res ]
        variances = [ x.variance for x in var_by_res ]
         
        # Get list of residues to keep under the different intevals
        truncation_levels = []
        truncation_variances = []
        truncation_residues = []
        truncation_residue_idxs = []
        for start in start_idxs:
            percent = int(round(float(start + 1) / float(length) * 100))
            residues = resseq_all[:start + 1]
            idxs = resseq_all[:start + 1]
            idxs = idxs_all[:start + 1]
            thresh = variances[start]  # For the threshold we take the threshold of the most variable residue
            truncation_variances.append(thresh)
            truncation_levels.append(percent)
            # print "GOT PERCENT,THRESH ",percent,thresh
            # print "residues ",residues
            truncation_residues.append(sorted(residues))
            truncation_residue_idxs.append(sorted(idxs))
                 
        return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs
    
    def _calculate_residues_thresh(self, var_by_res, percent_interval):
        """Txxx
        """

        # calculate the thresholds
        truncation_variances = self.generate_thresholds(var_by_res, percent_interval)

        # We run in reverse as that's how the original code worked
        truncation_residues = []
        truncation_residue_idxs = []
        truncation_levels = []
        lt = len(truncation_variances)
        for i, truncation_threshold in enumerate(truncation_variances):
            
            truncation_level = lt - i  # as going backwards
            truncation_levels.append(truncation_level)
            
            # Get a list of the indexes of the residues to keep
            to_keep = [ x.resSeq for x in var_by_res if x.variance <= truncation_threshold ]
            to_keep_idxs = [ x.idx for x in var_by_res if x.variance <= truncation_threshold ]
            truncation_residues.append(to_keep)
            truncation_residue_idxs.append(to_keep_idxs)
        
        # We went through in reverse so put things the right way around
        truncation_levels.reverse()
        truncation_variances.reverse()
        truncation_residues.reverse()
        truncation_residue_idxs.reverse()
        return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs
        
    def cluster_models(self,
                       models=None,
                       cluster_method=None,
                       num_clusters=None,
                       cluster_exe=None,
                       import_cluster=False,
                       cluster_dir=None,
                       nproc=1,
                       max_cluster_size=200
                       ):
        clusters = []
        clusters_data = []
        if import_cluster:
            if not os.path.isdir(cluster_dir): raise RuntimeError, "Import cluster cannot find directory: {0}".format(cluster_dir)
            cluster_models = glob.glob(os.path.join(cluster_dir, "*.pdb"))
            if not cluster_models: raise RuntimeError, "Import cluster cannot find pdbs in directory: {0}".format(cluster_dir)
            # Data on the cluster
            cluster_data = self.create_dict()
            cluster_data['cluster_num'] = 1
            cluster_data['cluster_centroid'] = cluster_models[0]
            cluster_data['cluster_num_models'] = len(cluster_models)
            cluster_data['cluster_method'] = "import"
            cluster_data['num_clusters'] = 1
            clusters_data.append(cluster_data)
            clusters.append(cluster_models)          
               
        elif cluster_method == "spicker":
            # Spicker Alternative for clustering
            self.logger.info('* Running SPICKER to cluster models *')
            spicker_rundir = os.path.join(self.work_dir, 'spicker')
            spickerer = spicker.Spickerer(spicker_exe=cluster_exe)
            spickerer.cluster(models, run_dir=spicker_rundir)
            self.logger.debug(spickerer.results_summary())
            
            ns_clusters=len(spickerer.results)
            if ns_clusters == 0: raise RuntimeError,"No clusters returned by SPICKER"
            if ns_clusters < num_clusters:
                self.logger.critical('Requested {0} clusters but SPICKER only found {0} so using {1} clusters'.format(num_clusters,ns_clusters))
                num_clusters=ns_clusters
            
            for i in range(num_clusters):
                # We truncate the list of models to max_cluster_size. This probably needs to be redone, because as the models are ordered by their similarity
                # to the cluster centroid, we automatically select the 200 most similar to the centroid. However if the cluster is large and the models similar
                # then when theseus calculates the variances, the variances will be representative of the 200, but might not show how the models vary throughout
                # the whole cluster, which could provide better information for truncating the models.
                cluster = spickerer.results[i].pdbs[0:max_cluster_size]
                clusters.append(cluster)
                # Data on the models
                cluster_data = self.create_dict()
                d = spickerer.results[i]
                cluster_data['cluster_num'] = i + 1
                cluster_data['cluster_centroid'] = d.cluster_centroid
                cluster_data['cluster_num_models'] = len(cluster)
                cluster_data['cluster_method'] = cluster_method
                cluster_data['num_clusters'] = num_clusters
                clusters_data.append(cluster_data)
        
        elif cluster_method == "fast_protein_cluster":
            fpc = fast_protein_cluster.FPC()
            SCORE_TYPE = 'rmsd'
            CLUSTER_METHOD = 'kmeans'
            self.logger.info('Running fast_protein_cluster with: score_type: {0} cluster_method: {1}'.format(SCORE_TYPE,
                                                                                                             CLUSTER_METHOD))
            fpc_rundir = os.path.join(self.work_dir, 'fast_protein_cluster')
            self.logger.info('fast_protein_cluster running in directory: {0}'.format(fpc_rundir))
            clusters, clusters_data = fpc.cluster(models=models,
                                                num_clusters=num_clusters,
                                                score_type=SCORE_TYPE,
                                                cluster_method=CLUSTER_METHOD,
                                                work_dir=fpc_rundir,
                                                fpc_exe=cluster_exe,
                                                nproc=nproc,
                                                max_cluster_size=max_cluster_size)
        else:
            raise RuntimeError, 'Unrecognised clustering method: {0}'.format(cluster_method)

        return clusters, clusters_data
    
    def create_dict(self):
        """Create an empty dictionary
        Not strictly necessary but it's a place to remember what we capture
        """
        d = {}
        d['cluster_method'] = None
        d['num_clusters'] = None
        d['cluster_num'] = None
        d['cluster_centroid'] = None
        d['cluster_num_models'] = None
        
        # truncation info
        d['truncation_level'] = None
        d['percent_truncation'] = None
        d['truncation_method'] = None
        d['truncation_residues'] = None
        d['truncation_dir'] = None
        d['truncation_variance'] = None
        d['num_residues'] = None

        # subclustering info
        d['subcluster_num_models'] = None
        d['subcluster_radius_threshold'] = None
        d['subcluster_centroid_model'] = None
    
        # ensemble info
        d['name'] = None
        d['side_chain_treatment'] = None
        d['ensemble_num_atoms'] = None
        d['ensemble_pdb'] = None  # path to the ensemble file
        
        return d
    
    def edit_side_chains(self, raw_ensemble, raw_ensemble_data, ensembles_directory, homologs=False, side_chain_treatments=SIDE_CHAIN_TREATMENTS):
        assert os.path.isdir(ensembles_directory), "Cannot find ensembles directory: {0}".format(ensembles_directory)
        ensembles = []
        ensembles_data = []
        if side_chain_treatments is None: side_chain_treatments=[UNMODIFIED]
        for sct in side_chain_treatments:
            ensemble_data = copy.copy(raw_ensemble_data)
            ensemble_data['side_chain_treatment'] = sct
            if homologs:
                ensemble_data['name'] = 'e{0}_{1}'.format(ensemble_data['truncation_level'], sct)
            else:
                ensemble_data['name'] = 'c{0}_t{1}_r{2}_{3}'.format(ensemble_data['cluster_num'],
                                                                   ensemble_data['truncation_level'],
                                                                   ensemble_data['subcluster_radius_threshold'],
                                                                   sct)
            # create filename based on name and side chain treatment
            # fpath = ample_util.filename_append(raw_ensemble,astr=sct, directory=ensembles_directory)
            fpath = os.path.join(ensembles_directory, "{0}.pdb".format(ensemble_data['name']))
            
            # Create the files
            if sct == ALLATOM or sct == UNMODIFIED:
                # For all atom just copy the file
                shutil.copy2(raw_ensemble, fpath)
            elif sct == RELIABLE:
                pdb_edit.reliable_sidechains(raw_ensemble, fpath)
            elif sct == POLYALA:
                pdb_edit.backbone(raw_ensemble, fpath)
            else:
                raise RuntimeError, "Unrecognised side_chain_treatment: {0}".format(sct)
            
            # Count the number of atoms in the ensemble-only required for benchmark mode
            natoms, nresidues = pdb_edit.num_atoms_and_residues(fpath, first=True)
            
            # Process ensemble data
            ensemble_data['ensemble_pdb'] = fpath
            ensemble_data['ensemble_num_atoms'] = natoms
            # check
            assert ensemble_data['num_residues'] == nresidues, "Unmatching number of residues: {0} : {1} \n{2}".format(ensemble_data['num_residues'],
                                                                                                                       nresidues,
                                                                                                                       raw_ensemble_data)
            ensembles.append(fpath)
            ensembles_data.append(ensemble_data)
                
        return ensembles, ensembles_data
  
    def generate_ensembles(self, models,
                           cluster_method=None,
                           cluster_exe=None,
                           num_clusters=None,
                           import_cluster=False,
                           cluster_dir=None,
                           percent_truncation=None,
                           truncation_method=None,
                           truncation_pruning=None,
                           ensembles_directory=None,
                           work_dir=None,
                           nproc=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS):
        
        # Work dir set each time
        if not work_dir:
            raise RuntimeError, "Need to set work_dir!"
        self.work_dir = work_dir
        
        if not cluster_method:
            cluster_method = self.cluster_method
        if not cluster_exe:
            cluster_exe = self.cluster_exe
        if not num_clusters:
            num_clusters = self.num_clusters
        if not percent_truncation:
            percent_truncation = self.percent_truncation
        if not truncation_method:
            truncation_method = self.truncation_method
        if not truncation_pruning:
            truncation_pruning = self.truncation_pruning
        if not ensembles_directory:
            self.ensembles_directory = os.path.join(work_dir, "ensembles")
        else:
            self.ensembles_directory = ensembles_directory
        
        if not import_cluster and not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError, "Problem reading models given to Ensembler: {0}".format(models) 
        
        self.logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory):
            os.mkdir(self.ensembles_directory)

        self.ensembles = []
        self.ensembles_data = []
        for cluster, cluster_data in zip(*self.cluster_models(models=models,
                                                              cluster_method=cluster_method,
                                                              num_clusters=num_clusters,
                                                              cluster_exe=cluster_exe,
                                                              import_cluster=import_cluster,
                                                              cluster_dir=cluster_dir,
                                                              nproc=nproc)):
            if len(cluster) < 2:
                self.logger.info("Cannot truncate cluster {0} as < 2 models!".format(cluster_data['cluster_num']))
                continue
            self.logger.info('Processing cluster: {0}'.format(cluster_data['cluster_num']))
            for truncated_models, truncated_models_data in zip(*self.truncate_models(cluster,
                                                                                     cluster_data,
                                                                                     truncation_method=truncation_method,
                                                                                     truncation_pruning=truncation_pruning,
                                                                                     percent_truncation=percent_truncation)):
                for subcluster, subcluster_data in zip(*self.subcluster_models(truncated_models,
                                                                               truncated_models_data,
                                                                               subcluster_program=self.subcluster_program,
                                                                               subcluster_exe=self.subcluster_program,
                                                                               ensemble_max_models=self.ensemble_max_models)):
                    for ensemble, ensemble_data in zip(*self.edit_side_chains(subcluster, subcluster_data, self.ensembles_directory)):
                        self.ensembles.append(ensemble)
                        self.ensembles_data.append(ensemble_data)
        
        return self.ensembles

    def generate_ensembles_homologs(self,
                                    models,
                                    alignment_file=None,
                                    percent_truncation=None,
                                    truncation_method=None,
                                    ensembles_directory=None,
                                    work_dir=None,
                                    nproc=None,
                                    homolog_aligner=None,
                                    mustang_exe=None,
                                    gesamt_exe=None,
                                    side_chain_treatments=SIDE_CHAIN_TREATMENTS):
        
        # Work dir set each time
        if not work_dir: raise RuntimeError, "Need to set work_dir!"
        self.work_dir = work_dir
        
        if not percent_truncation:
            percent_truncation = self.percent_truncation
        if not truncation_method:
            truncation_method = self.truncation_method
        if not ensembles_directory:
            self.ensembles_directory = os.path.join(work_dir, "ensembles")
        else:
            self.ensembles_directory = ensembles_directory
        
        if not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError, "Problem reading models given to Ensembler: {0}".format(models) 
        
        self.logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)
        
        # standardise all the models
        std_models_dir = os.path.join(work_dir, "std_models")
        os.mkdir(std_models_dir)
        std_models = []
        for m in models:
            std_model = ample_util.filename_append(m, 'std', std_models_dir)
            pdb_edit.standardise(pdbin=m, pdbout=std_model, del_hetatm=True)
            std_models.append(std_model)
        
        if not alignment_file:
            if homolog_aligner == 'mustang':
                self.logger.info("Generating alignment file with mustang_exe: {0}".format(mustang_exe))
                alignment_file = align_mustang(std_models, mustang_exe=mustang_exe, work_dir=self.work_dir)
            elif homolog_aligner == 'gesamt':
                self.logger.info("Generating alignment file with gesamt_exe: {0}".format(gesamt_exe))
                alignment_file = align_gesamt(std_models, gesamt_exe=gesamt_exe, work_dir=self.work_dir)
            else:
                raise RuntimeError, "Unknown homolog_aligner: {0}".format(homolog_aligner)
            self.logger.info("Generated alignment file: {0}".format(alignment_file))
        else:
            self.logger.info("Using alignment file: {0}".format(alignment_file))
            
        # Now truncate and create ensembles - as standard ample, but with no subclustering
        self.ensembles = []
        self.ensembles_data = []
        for truncated_models, truncated_models_data in zip(*self.truncate_models(std_models,
                                                                                 truncation_method=truncation_method,
                                                                                 truncation_pruning=None,
                                                                                 percent_truncation=percent_truncation,
                                                                                 homologs=True,
                                                                                 alignment_file=alignment_file
                                                                                 )):
            tlevel = truncated_models_data['truncation_level']
            ensemble_dir = os.path.join(truncated_models_data['truncation_dir'],
                                        "ensemble_{0}".format(tlevel))
            os.mkdir(ensemble_dir)
            os.chdir(ensemble_dir)
             
            # Need to create an alignment file for theseus
            basename = "e{0}".format(tlevel)
            pre_ensemble = self.align_models(truncated_models, basename=basename, work_dir=ensemble_dir, homologs=True)
            if not pre_ensemble:
                self.logger.critical("Skipping ensemble {0} due to error with Theseus".format(basename))
                continue
            pre_ensemble_data = copy.copy(truncated_models_data)
             
            for ensemble, ensemble_data in zip(*self.edit_side_chains(pre_ensemble,
                                                                      pre_ensemble_data,
                                                                      self.ensembles_directory,
                                                                      homologs=True,
                                                                      side_chain_treatments=side_chain_treatments)):
                self.ensembles.append(ensemble)
                self.ensembles_data.append(ensemble_data)
        
        return self.ensembles

    def generate_thresholds(self, var_by_res, percent_interval):
        """
        This is the original method developed by Jaclyn and used in all work until November 2014 (including the coiled-coil paper)
        
        Calculate the residue variance thresholds that will keep self.percent_interval residues for each truncation level
        """
        #--------------------------------
        # choose threshold type
        #-------------------------------
        FIXED_INTERVALS = False
        if FIXED_INTERVALS:
            self.thresholds = [ 1, 1.5, 2 , 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8 ]
            self.logger.debug("Got {0} thresholds: {1}".format(len(self.thresholds), self.thresholds))
            return

        # List of variances ordered by residue index
        var_list = [ x.variance for x in var_by_res]
        length = len(var_list)
        if length == 0:
            msg = "Error generating thresholds, got len: {0}".format(length)
            self.logger.critical(msg)
            raise RuntimeError, msg

        # How many residues should fit in each bin
        # NB - Should round up not down with int!
        chunk_size = int((float(length) / 100) * float(percent_interval))
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            self.logger.critical(msg)
            raise RuntimeError, msg

        # # try to find intervals for truncation
        truncation_thresholds = self._generate_thresholds(var_list, chunk_size)
        
        # Jens' new untested method
        # truncation_thresholds=self._generate_thresholds2(var_list, chunk_size)
        
        self.logger.debug("Got {0} thresholds: {1}".format(len(truncation_thresholds), truncation_thresholds))

        return truncation_thresholds

    def _generate_thresholds(self, values, chunk_size):
        """Jaclyn's threshold method
        """
        try_list = copy.deepcopy(values)
        try_list.sort()
        # print "try_list ",try_list

        # print list(chunks(try_list, int(chunk_size)))
        # For chunking list
        def chunks(a_list, chunk_size):
            for i in xrange(0, len(a_list), chunk_size):
                yield a_list[i:i + chunk_size ]

        thresholds = []
        for x in list(chunks(try_list, chunk_size)):
            # print x, x[-1]
            # For some cases, multiple residues share the same variance so we don't create a separate thereshold
            if x[-1] not in thresholds:
                thresholds.append(x[-1])
                
        return thresholds

    def _generate_thresholds2(self, values, chunk_size):
        """
        This is Jens's update to Jaclyn's method that groups the residues by variances so that we split
        them by variance, and try and fit chunk_size in each bin. Previously we tried to split by variance but didn't
        group the residues by variance, so the same variance bin could cover multiple residue groups. 
        """

        # Create tuple mapping values to counts
        data = [(i, values.count(i)) for i in sorted(set(values), reverse=True)]

        thresholds = []
        counts = []
        first = True
        for variance, count in data:
            if first or counts[-1] + count > chunk_size:
                thresholds.append(variance)
                counts.append(count)
                if first: first = False
            else:
                # thresholds[-1]=variance
                counts[-1] += count

        thresholds.sort()
        return thresholds
    
    def prune_residues(self, residues, chunk_size=1, allowed_gap=2):
        """Remove any residues that are < chunk_size where the gap before and after is > allowed_gap"""
        
        assert chunk_size > 0 and allowed_gap > 0, "chunk_size and allowed_gap must be > 0!: {0} {1}".format(chunk_size, allowed_gap)
        
        if not len(residues): return residues, None
        lenr = len(residues)
        if lenr <= chunk_size:
            return [], residues
        
        # Build up a list of residues to remove
        to_remove = []
        start = residues[0]
        last = residues[0]
        this_residue = None
        last_chunk_end = residues[0] - (allowed_gap + 1)  # make sure starting gap is bigger than allowed 
        i = 1
        idxLast = lenr - 1
        while i <= idxLast:
            this_residue = residues[i]
            if i == idxLast or this_residue != last + 1:
                if i == idxLast:  # Need to fiddle things at the end
                    if this_residue != last + 1:  # handle last if a singleton
                        start = this_residue
                        last_chunk_end = last
                    last = this_residue
                    postgap = allowed_gap + 1  # at end so just larger then gap
                else:
                    postgap = (this_residue - last) - 1
                
                pregap = (start - last_chunk_end) - 1
                this_chunk_size = (last - start) + 1
                
                # remove if it satisfies the requirements
                if (this_chunk_size <= chunk_size and pregap >= allowed_gap and postgap >= allowed_gap):
                    chunk = [x for x in range(start, last + 1)]
                    to_remove += chunk
                
                # reset start and last_chunk_end
                start = this_residue
                last_chunk_end = last
                
            # update loop
            last = this_residue
            i += 1
        
        # Remove the chunks and return
        if len(to_remove):
            return [r for r in residues if r not in to_remove], to_remove
        else:
            return residues, None

    def subcluster_models(self,
                           truncated_models,
                           truncated_models_data,
                           subcluster_program=None,
                           subcluster_exe=None,
                           ensemble_max_models=None):
        
        if self.subcluster_method == "ORIGINAL":
            f = self.subcluster_models_fixed_radii
        elif self.subcluster_method == "FLOATING_RADII":
            f = self.subcluster_models_floating_radii
        else:
            assert False
        return f(truncated_models,
                 truncated_models_data,
                 subcluster_program,
                 subcluster_exe,
                 ensemble_max_models)
        

            
    def subcluster_models_fixed_radii(self,
                                      truncated_models,
                                      truncated_models_data,
                                      subcluster_program=None,
                                      subcluster_exe=None,
                                      ensemble_max_models=None,
                                      radius_thresholds=None
                                      ):
        
        # Theseus only works with > 3 residues
        if truncated_models_data['num_residues'] <= 2: return [], []
        
        if not radius_thresholds: radius_thresholds = self.subcluster_radius_thresholds
        ensembles = []
        ensembles_data = []
        
        # Use first model to get data on level
        cluster_num = truncated_models_data['cluster_num']
        truncation_level = truncated_models_data['truncation_level']
        truncation_dir = truncated_models_data['truncation_dir']
        
        # Make sure everyting happens in the truncation directory
        owd = os.getcwd()
        os.chdir(truncation_dir)
            
        # Run maxcluster to generate the distance matrix
        if subcluster_program == 'maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else:
            assert False
        clusterer.generate_distance_matrix(truncated_models)
        # clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging

        # Loop through the radius thresholds
        previous_clusters = []
        for radius in radius_thresholds:
            self.logger.debug("subclustering models under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            if not cluster_files:
                self.logger.debug("Skipping radius {0} as no files clustered in directory {1}".format(radius, truncation_dir))
                continue
                
            self.logger.debug("Clustered {0} files".format(len(cluster_files)))
            cluster_files = self._slice_subcluster(cluster_files, previous_clusters, ensemble_max_models, radius, radius_thresholds)
            if not cluster_files:
                self.logger.debug('Could not create different cluster for radius {0} in directory: {1}'.format(radius, truncation_dir))
                continue
            
            # Remember this cluster so we don't create duplicate clusters
            previous_clusters.append(cluster_files)

            # Got files so create the directories
            subcluster_dir = os.path.join(truncation_dir, 'subcluster_{0}'.format(radius))
            os.mkdir(subcluster_dir)
            os.chdir(subcluster_dir)
            basename = 'c{0}_t{1}_r{2}'.format(cluster_num, truncation_level, radius)
            
            # List of files for reference
            with open(os.path.join(subcluster_dir, "{0}.list".format(basename)), 'w') as f:
                for m in cluster_files: f.write(m + "\n")
                f.write("\n")
            
            cluster_file = self.align_models(cluster_files, work_dir=subcluster_dir)
            if not cluster_file:
                msg = "Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename, subcluster_dir)
                self.logger.critical(msg)
                continue
             
            ensemble = os.path.join(subcluster_dir, basename + '.pdb')
            shutil.move(cluster_file, ensemble)

            # The data we've collected is the same for all pdbs in this level so just keep using the first  
            ensemble_data = copy.copy(truncated_models_data)
            ensemble_data['subcluster_num_models'] = len(cluster_files)
            ensemble_data['subcluster_radius_threshold'] = radius
            ensemble_data['ensemble_pdb'] = ensemble

            # Get the centroid model name from the list of files given to theseus - we can't parse
            # the pdb file as theseus truncates the filename
            ensemble_data['subcluster_centroid_model'] = os.path.abspath(cluster_files[0])
            
            ensembles.append(ensemble)
            ensembles_data.append(ensemble_data)
        
        # back to where we started
        os.chdir(owd)
        
        return ensembles, ensembles_data

    def _slice_subcluster(self, cluster_files, previous_clusters, ensemble_max_models, radius, radius_thresholds):
        """Select a unique set of models from a subcluster of models.
        """
        len_cluster = len(cluster_files)
        if not len_cluster: return None
        len_radius_thresholds = len(radius_thresholds)
        if len_cluster <= ensemble_max_models:
            if cluster_files not in previous_clusters: return cluster_files
            else: return None
        
        if len_cluster > ensemble_max_models:
            idx = radius_thresholds.index(radius)
            selected = cluster_files[:ensemble_max_models]
            if idx == 0 or selected not in previous_clusters: return selected
            
            # Here we have more models then we need, but the first slice has already been selected
            # we therefore need to select another slice
            
            # If last radius threshold, just take the slice to the end
            if idx + 1 == len_radius_thresholds:
                start = len_cluster - ensemble_max_models
                selected = cluster_files[start:]
                if selected not in previous_clusters:
                    return selected
                else:
                    return None
            
            # Work out how many residues are extra
            remainder = len_cluster - ensemble_max_models
            
            # Use the position of the radius in the list of radii to work out where to start this slice
            prop = float(idx) / float(len(radius_thresholds) - 1)  # -1 as the first is always at the start
            
            # Work out how many residues in to the remainder to start
            start = int(round(float(remainder) * prop))
            selected = cluster_files[start :  start + ensemble_max_models]
            if selected and selected not in previous_clusters:
                    return selected
            else:
                return None
        
        return None

    def subcluster_models_floating_radii(self,
                                         truncated_models,
                                         truncated_models_data,
                                         subcluster_program=None,
                                         subcluster_exe=None,
                                         ensemble_max_models=None):
        self.logger.info("subclustering with floaing radii")

        # Run maxcluster to generate the distance matrix
        if subcluster_program == 'maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else: assert False
        clusterer.generate_distance_matrix(truncated_models)
        # clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging
        
        subclusters = []
        subclusters_data = []
        clusters = []
        radii = []
        len_truncated_models = len(truncated_models)
        for i in range(len(self.subcluster_radius_thresholds)):
            radius = None
            nmodels = None
            if i > 0 and radii[i - 1] > self.subcluster_radius_thresholds[i]:
                radius = radii[i - 1]
                nmodels = len(clusters[i - 1])
                cluster_files, radius = self._subcluster_nmodels(nmodels, radius, clusterer, direction='up', increment=1)
            else:
                radius = self.subcluster_radius_thresholds[i]
                cluster_files = clusterer.cluster_by_radius(radius)
            
            if cluster_files:
                cluster_files = tuple(sorted(cluster_files))  # Need to sort so that we can check if we've had this lot before
                cluster_size = len(cluster_files)
            else:
                cluster_files = []
                cluster_size = 0

            if radius in radii or cluster_size == 0:
                # Increase radius till we have one more than the last one
                if cluster_size == 0:
                    nmodels = 2
                else:
                    radius = radii[i - 1]
                    nmodels = len(clusters[i - 1]) + 1
                cluster_files, radius = self._subcluster_nmodels(nmodels, radius, clusterer, direction='up', increment=1)
                cluster_files = sorted(cluster_files)
            elif cluster_size >= ensemble_max_models or cluster_files in clusters:
                # Randomly pick ensemble_max_models
                cluster_files = self._pick_nmodels(cluster_files, clusters, ensemble_max_models)
                if not cluster_files:
                    self.logger.debug('Could not cluster files under radius: {0} - could not find different models'.format(radius))
                    break
            
            # Need to check in case we couldn't cluster under this radius
            if cluster_size == 0 or radius in radii:
                self.logger.debug('Could not cluster files under radius: {0} - got {1} files'.format(radius, len(cluster_files)))
                break
            
            self.logger.debug('Subclustering {0} files under radius {1}'.format(cluster_size, radius))
            try:
                cluster_ensemble, data = self._subcluster_radius(list(cluster_files), radius, truncated_models_data)
            except RuntimeError:
                self.logger.debug('Could not cluster files under radius: {0} as theseus failed'.format(radius, len(cluster_files)))
                # If theseus fails, we just move
                break
            
            subclusters.append(cluster_ensemble)
            subclusters_data.append(data)
            clusters.append(tuple(cluster_files))  # append as tuple so it is hashable
            radii.append(radius)
            if cluster_size == len_truncated_models: break
            
        return subclusters, subclusters_data

    def _pick_nmodels(self, models, clusters, ensemble_max_models):
        MAXTRIES = 50
        tries = 0
        clusters = set(clusters)
        nmodels = min(len(models), ensemble_max_models)
        while True:
            subcluster = random.sample(models, nmodels)
            subcluster = tuple(sorted(subcluster))
            if subcluster not in clusters: break
            tries += 1
            if tries >= MAXTRIES: return None
        return subcluster
    
    def _subcluster_nmodels(self, nmodels, radius, clusterer, direction, increment):
        """
        """
        MINRADIUS = 0.0001
        MAXRADIUS = 100
        subcluster_models = clusterer.cluster_by_radius(radius)
        if subcluster_models:
            len_models = len(subcluster_models)
        else:
            len_models = 0
        self.logger.debug("_subcluster_nmodels: {0} {1} {2} {3} {4}".format(len_models, nmodels, radius, direction, increment))
        if len_models == nmodels or radius < MINRADIUS or radius > MAXRADIUS:
            self.logger.debug("_subcluster_nmodels returning: nmodels: {0} radius: {1}".format(len_models, radius))
            return subcluster_models, radius
        
        def lower_increment(increment):
            increment = increment / float(10)
            if increment <= 0.00001: raise RuntimeError, "increment out of bounds"
            return increment
        
        # Am sure the logic could be improved here, but it seems to  work
        try:
            if len_models > nmodels:
                # If we have more models than we want and we are increasing the radius, we've overshot, so we need to
                # decrease the radius but by a smaller increment
                # If the radius is the same as the increment, we need to decrease the incrememnt before we subtract it
                # as both of the above require decreasing the increment we have one test and just change the direction
                # for the overshoot
                if direction == 'up' or abs(radius - increment) < 0.0000001:
                    if direction == 'up': direction = 'down'
                    increment = lower_increment(increment)
                radius -= increment
            elif len_models < nmodels:
                if direction == 'down' :
                    direction = 'up'
                    increment = lower_increment(increment)
                radius += increment
        except RuntimeError:
            # Can't get a match so just return what we have
            self.logger.debug("_subcluster_nmodels exceeded increment. Returning: nmodels: {0} radius: {1}".format(len(subcluster_models), radius))
            return subcluster_models, radius
            
        return self._subcluster_nmodels(nmodels, radius, clusterer, direction, increment)
    
    def _subcluster_radius(self, models, radius, truncated_models_data):
        # Extract data from dictionary
        cluster_num = truncated_models_data['cluster_num']
        truncation_level = truncated_models_data['truncation_level']
        truncation_dir = truncated_models_data['truncation_dir']

        # Got files so create the directories
        subcluster_dir = os.path.join(truncation_dir, 'subcluster_{0}'.format(radius))
        os.mkdir(subcluster_dir)
        os.chdir(subcluster_dir)

        basename = 'c{0}_t{1}_r{2}'.format(cluster_num, truncation_level, radius)
        cluster_file = self.align_models(models)
        if not cluster_file:
            msg = "Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename,
                                                                                                subcluster_dir)
            self.logger.critical(msg)
            raise RuntimeError, msg
        
        ensemble = os.path.join(subcluster_dir, basename + '.pdb')
        shutil.move(cluster_file, ensemble)

        # The data we've collected is the same for all pdbs in this level so just keep using the first  
        subcluster_data = copy.copy(truncated_models_data)
        subcluster_data['subcluster_num_models'] = len(models)
        subcluster_data['subcluster_radius_threshold'] = radius
        subcluster_data['ensemble_pdb'] = ensemble

        # Get the centroid model name from the list of files given to theseus - we can't parse
        # the pdb file as theseus truncates the filename
        subcluster_data['subcluster_centroid_model'] = os.path.abspath(models[0])
        return ensemble, subcluster_data
  
    def truncate_models(self,
                        models,
                        models_data={},
                        truncation_method=None,
                        percent_truncation=None,
                        truncation_pruning='none',
                        homologs=False,
                        alignment_file=None
                        ):
        
        assert len(models) > 1, "Cannot truncate as < 2 models!"
        assert truncation_method and percent_truncation, "Missing arguments: {0}".format(truncation_method)

        # Create the directories we'll be working in
        if homologs: truncate_dir = os.path.join(self.work_dir, 'truncate')
        else: truncate_dir = os.path.join(self.work_dir, 'truncate_{0}'.format(models_data['cluster_num']))
        os.mkdir(truncate_dir)
        os.chdir(truncate_dir)
        
        # Calculate variances between pdb - and align them if necessary
        run_theseus = theseus.Theseus(work_dir=truncate_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.align_models(models, homologs=homologs, alignment_file=alignment_file)
        except RuntimeError,e:
            self.logger.critical(e)
            return [],[]
            
        var_by_res = run_theseus.var_by_res()
        if not len(var_by_res) > 0:
            msg = "Error reading residue variances!"
            self.logger.critical(msg)
            raise RuntimeError, msg
        
        # Need to trim the aligned models down to core
        if homologs: models = model_core_from_alignment(run_theseus.aligned_models, alignment_file)
            
        self.logger.info('Using truncation method: {0}'.format(truncation_method))
        # Calculate which residues to keep under the different methods
        truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = None, None, None, None
        if truncation_method == 'percent':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = self._calculate_residues_percent(var_by_res, percent_truncation)
        elif truncation_method == 'thresh':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = self._calculate_residues_thresh(var_by_res, percent_truncation)
        elif truncation_method == 'focussed':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = self._calculate_residues_focussed(var_by_res)
        else:
            raise RuntimeError, "Unrecognised ensembling mode: {0}".format(truncation_method)
        
        self.truncation_levels = truncation_levels # save so we can put in results dict
        self.truncation_variances = truncation_variances # save so we can put in results dict
        self.truncation_nresidues = [len(r) for r in truncation_residues] # save so we can put in results dict

        truncated_models = []
        truncated_models_data = []
        pruned_residues = None
        for tlevel, tvar, tresidues, tresidue_idxs in zip(truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs):
            # Prune singletone/doubletone etc. residues if required
# This code is currently disabled as it relies on resseqs and we now use residue indexes to allow us to work with pdbs with different reseqs
#             self.logger.debug("truncation_pruning: {0}".format(truncation_pruning))
#             if truncation_pruning=='single':
#                 tresidues,pruned_residues=self.prune_residues(tresidues, chunk_size=1, allowed_gap=2)
#                 if pruned_residues: self.logger.debug("prune_residues removing: {0}".format(pruned_residues))
#             elif truncation_pruning=='none':
#                 pass
#             else:
#                 raise RuntimeError,"Unrecognised truncation_pruning: {0}".format(truncation_pruning)
            
            # Skip if there are no residues
            if not tresidue_idxs:
                self.logger.debug("Skipping truncation level {0} with variance {1} as no residues".format(tlevel, tvar))
                continue
            
            trunc_dir = os.path.join(truncate_dir, 'tlevel_{0}'.format(tlevel))
            os.mkdir(trunc_dir)
            self.logger.info('Truncating at: {0} in directory {1}'.format(tlevel, trunc_dir))

            # list of models for this truncation level
            level_models = []
            for infile in models:
                pdbname = os.path.basename(infile)
                pdbout = os.path.join(trunc_dir, pdbname)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdb_edit.select_residues(pdbin=infile, pdbout=pdbout, tokeep_idx=tresidue_idxs)
                level_models.append(pdbout)
            
            # Add the model
            truncated_models.append(level_models)

            # Add the data
            model_data = copy.copy(models_data)
            model_data['truncation_level'] = tlevel
            model_data['truncation_variance'] = tvar
            model_data['truncation_residues'] = tresidues
            model_data['num_residues'] = len(tresidues)
            model_data['truncation_dir'] = trunc_dir
            model_data['percent_truncation'] = percent_truncation
            model_data['truncation_method'] = truncation_method
            model_data['truncation_pruning'] = truncation_pruning
            model_data['pruned_residues'] = pruned_residues
            
            truncated_models_data.append(model_data)
            
        return truncated_models, truncated_models_data

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


        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return

    def testPruneResidues(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()

        ensembler = Ensembler()
        
        residues = []
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, None)
          
        residues = [1]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, [1])
         
        residues = [1, 2]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, [1, 2])
         
        residues = [1, 2, 4]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 4])
        self.assertEqual(pruned, None)
         
        # Big enough gap 
        residues = [1, 2, 3, 4, 8, 13, 14, 15]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 3, 4, 13, 14, 15])
        self.assertEqual(pruned, [8])
           
        residues = [1, 2, 3, 4, 8, 9, 13, 14, 15]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 3, 4, 13, 14, 15])
        self.assertEqual(pruned, [8, 9])
           
        residues = [1, 2, 3, 4, 8, 9, 10, 13, 14, 15]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, residues)
        self.assertEqual(pruned, None)
           
        # end gap not big enough
        residues = [1, 2, 3, 4, 8, 9, 11, 12, 13]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, residues)
        self.assertEqual(pruned, None)
           
        # Lone residue at start
        residues = [1, 11, 12, 13]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [11, 12, 13])
        self.assertEqual(pruned, [1])
        
        # Lone residue at end
        residues = [11, 12, 13, 19]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [11, 12, 13])
        self.assertEqual(pruned, [19])
          
        # Mixed
        residues = [1, 3, 4, 7, 10, 11, 13, 15, 16, 19]
        pres, pruned = ensembler.prune_residues(residues, chunk_size=1, allowed_gap=2)
        self.assertEqual(pres, [1, 3, 4, 10, 11, 13, 15, 16])
        self.assertEqual(pruned, [7, 19])

        return

    def testThresholds(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()

        ensembler = Ensembler()
        
        ensembler.work_dir = os.path.join(self.tests_dir, "genthresh1")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        
        percent_interval = 5
        mdir = os.path.join(self.testfiles_dir, "models")
        cluster_models = glob.glob(mdir + os.sep + "*.pdb")
        run_theseus = theseus.Theseus(theseus_exe=self.theseus_exe)
        run_theseus.align_models(cluster_models)
        var_by_res = run_theseus.var_by_res()
        thresholds = ensembler.generate_thresholds(var_by_res, percent_interval)
        
        self.assertEqual(30, len(thresholds), thresholds)
        reft = [0.02235, 0.041896, 0.08155, 0.085867, 0.103039, 0.185265, 0.7058, 1.772884, 3.152793,
              4.900255, 6.206563, 10.250043, 10.772177, 16.071345, 18.292366, 22.432564, 23.265938,
              25.348501, 27.37951, 35.877391, 37.227859, 41.759805, 49.037625, 50.992344, 56.091738,
              58.09387, 59.917999, 70.052007, 80.235358, 86.028994]
        
        self.assertTrue(all([ abs(r - c) < 0.0001 for r, c in zip(reft, thresholds)]),
                         "Wrong truncation thresholds: {0}".format(thresholds))
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testResiduesThresh(self):
        """Test we can calculate the original list of residues"""
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()
        
        # This test for percent
        ensembler.work_dir = os.path.join(self.tests_dir, "genthresh2")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        os.chdir(ensembler.work_dir)
        
        percent_interval = 5
        mdir = os.path.join(self.testfiles_dir, "models")
        cluster_models = glob.glob(mdir + os.sep + "*.pdb")
        run_theseus = theseus.Theseus(theseus_exe=self.theseus_exe)
        run_theseus.align_models(cluster_models)
        var_by_res = run_theseus.var_by_res()
        truncation_levels, truncation_variances, truncation_residues, truncation_residues_idxs = ensembler._calculate_residues_thresh(var_by_res, percent_interval)
        
        self.assertEqual(truncation_levels,
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        
        
        refv = [86.028994, 80.235358, 70.052007, 59.917999, 58.09387, 56.091738, 50.992344, 49.037625, 41.759805, 37.227859, 35.877391,
                27.37951, 25.348501, 23.265938, 22.432564, 18.292366, 16.071345, 10.772177, 10.250043, 6.206563, 4.900255, 3.152793,
                1.772884, 0.7058, 0.185265, 0.103039, 0.085867, 0.08155, 0.041896, 0.02235]
        self.assertTrue(all([ abs(r - c) < 0.0001 for r, c in zip(refv, truncation_variances)]),
                         "Wrong truncation variances: {0}".format(truncation_variances))

        residues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 56, 57],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 57],
                  [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
                  [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 53],
                  [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
                  [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
                  [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
                  [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45],
                  [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                  [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                  [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                  [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                  [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                  [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                  [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                  [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                  [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                  [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                  [25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                  [26, 27, 28, 30, 31, 32, 33, 34],
                  [26, 27, 30, 31, 32, 33],
                  [26, 27, 31, 32],
                  [31, 32]]
        self.assertEqual(residues, truncation_residues)

        shutil.rmtree(ensembler.work_dir)
        return
    
    def testResiduesFocussed(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()
        
        TheseusVariances = collections.namedtuple('TheseusVariances', ['idx', 'resName', 'resSeq', 'variance', 'stdDev', 'rmsd', 'core'])
        l = 160
        var_by_res = [ TheseusVariances(idx=i, resName='', resSeq=i, variance=float(i+1), stdDev=None, rmsd=None, core=None) for i in range(l) ]
        truncation_levels, truncation_variances, truncation_residues, truncation_residues_idxs = ensembler._calculate_residues_focussed(var_by_res)
         
        self.assertEqual(truncation_levels, [100, 93, 85, 78, 70, 63, 55, 48, 40, 33, 25, 23, 20, 18, 15, 13, 10, 8, 5, 3])
        self.assertEqual(truncation_variances, [160.0, 148.0, 136.0, 124.0, 112.0, 100.0, 88.0, 76.0, 64.0, 52.0, 40.0, 36.0, 32.0, 28.0, 24.0, 20.0, 16.0, 12.0, 8.0, 4.0])
        self.assertEqual(truncation_residues, truncation_residues_idxs)
        self.assertEqual(truncation_residues_idxs[-1], [0, 1, 2, 3])
       
        return
    
    def testResiduesPercent(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        ensembler.work_dir = os.path.join(self.tests_dir, "genthresh3")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        os.chdir(ensembler.work_dir)
        ensembler.percent_interval = 5
        mdir = os.path.join(self.testfiles_dir, "models")
        cluster_models = glob.glob(mdir + os.sep + "*.pdb")
        run_theseus = theseus.Theseus(theseus_exe=self.theseus_exe)
        run_theseus.align_models(cluster_models)
        var_by_res = run_theseus.var_by_res()
        truncation_levels, truncation_variances, truncation_residues, truncation_residues_idxs = ensembler._calculate_residues_percent(var_by_res, percent_interval=5)

        self.assertEqual(truncation_levels,
                         [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8])
        
        refv = [83.279358, 67.986091, 57.085407, 54.341361, 47.73422, 40.413985, 35.141765, 26.06671, 22.81897,
              21.187131, 14.889563, 9.55953, 6.484837, 4.263945, 1.41278, 0.504414, 0.204918, 0.135812, 0.081846]
        self.assertTrue(all([ abs(r - c) < 0.0001 for r, c in zip(refv, truncation_variances)]),
                         "Wrong truncation variances: {0}".format(truncation_variances))
        
        residues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 57],
                    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 53, 57],
                    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
                    [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
                    [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
                    [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                    [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                    [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                    [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                    [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                    [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],
                    [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
                    [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                    [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
                    [23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                    [26, 27, 28, 30, 31, 32, 33, 34],
                    [26, 27, 30, 31, 32]]

        for i in range(len(residues)):
            self.assertEqual(residues[i], truncation_residues[i], "Mismatching residues for level {0}\n{1}\n{2}".format(i, residues[i], truncation_residues[i]))

        shutil.rmtree(ensembler.work_dir)
        return

    def testClustering(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        ensembler.work_dir = os.path.join(self.tests_dir, "genthresh4")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")
        
        cluster_models, cluster_data = ensembler.cluster_models(models=models,
                                                             cluster_method='spicker',
                                                             num_clusters=3,
                                                             cluster_exe=self.spicker_exe)
        
        
        # This with spicker from ccp4 6.5.010 on osx 10.9.5
        # These should match with the test results from python/spicker.py
        names = sorted([os.path.basename(m) for m in cluster_models[0]])
        self.assertEqual(names,
                         sorted(['5_S_00000005.pdb', '4_S_00000005.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb',
                          '4_S_00000003.pdb', '3_S_00000006.pdb', '3_S_00000004.pdb', '2_S_00000005.pdb',
                          '2_S_00000001.pdb', '3_S_00000003.pdb', '1_S_00000005.pdb', '1_S_00000002.pdb',
                          '1_S_00000004.pdb']))
        
        d = cluster_data[2]
        self.assertEqual(os.path.basename(d['cluster_centroid']), '1_S_00000001.pdb')
        self.assertEqual(d['cluster_num_models'], 1)
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testEnsemblingPercent(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh5")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")

        num_clusters = 1
        cluster_method = 'spicker'
        percent_truncation = 5
        truncation_method = "percent"
        ensembles = ensembler.generate_ensembles(models,
                                                 cluster_method=cluster_method,
                                                 cluster_exe=self.spicker_exe,
                                                 num_clusters=num_clusters,
                                                 percent_truncation=percent_truncation,
                                                 truncation_method=truncation_method,
                                                 work_dir=work_dir)
        
        # Below tested with ccp4 6.5.010 on osx 10.9.5
        eref = ['c1_t100_r2_allatom.pdb', 'c1_t100_r2_polyAla.pdb', 'c1_t100_r2_reliable.pdb', 'c1_t100_r3_allatom.pdb',
                 'c1_t100_r3_polyAla.pdb', 'c1_t100_r3_reliable.pdb', 'c1_t19_r1_allatom.pdb', 'c1_t19_r1_polyAla.pdb',
                 'c1_t19_r1_reliable.pdb', 'c1_t24_r1_allatom.pdb', 'c1_t24_r1_polyAla.pdb', 'c1_t24_r1_reliable.pdb',
                 'c1_t29_r1_allatom.pdb', 'c1_t29_r1_polyAla.pdb', 'c1_t29_r1_reliable.pdb', 'c1_t34_r1_allatom.pdb',
                 'c1_t34_r1_polyAla.pdb', 'c1_t34_r1_reliable.pdb', 'c1_t39_r1_allatom.pdb', 'c1_t39_r1_polyAla.pdb',
                 'c1_t39_r1_reliable.pdb', 'c1_t44_r1_allatom.pdb', 'c1_t44_r1_polyAla.pdb', 'c1_t44_r1_reliable.pdb',
                 'c1_t44_r2_allatom.pdb', 'c1_t44_r2_polyAla.pdb', 'c1_t44_r2_reliable.pdb', 'c1_t49_r1_allatom.pdb',
                 'c1_t49_r1_polyAla.pdb', 'c1_t49_r1_reliable.pdb', 'c1_t49_r2_allatom.pdb', 'c1_t49_r2_polyAla.pdb',
                 'c1_t49_r2_reliable.pdb', 'c1_t54_r1_allatom.pdb', 'c1_t54_r1_polyAla.pdb', 'c1_t54_r1_reliable.pdb',
                 'c1_t54_r2_allatom.pdb', 'c1_t54_r2_polyAla.pdb', 'c1_t54_r2_reliable.pdb', 'c1_t59_r1_allatom.pdb',
                 'c1_t59_r1_polyAla.pdb', 'c1_t59_r1_reliable.pdb', 'c1_t59_r2_allatom.pdb', 'c1_t59_r2_polyAla.pdb',
                 'c1_t59_r2_reliable.pdb', 'c1_t64_r1_allatom.pdb', 'c1_t64_r1_polyAla.pdb', 'c1_t64_r1_reliable.pdb',
                 'c1_t64_r2_allatom.pdb', 'c1_t64_r2_polyAla.pdb', 'c1_t64_r2_reliable.pdb', 'c1_t69_r1_allatom.pdb',
                 'c1_t69_r1_polyAla.pdb', 'c1_t69_r1_reliable.pdb', 'c1_t69_r2_allatom.pdb', 'c1_t69_r2_polyAla.pdb',
                 'c1_t69_r2_reliable.pdb', 'c1_t75_r1_allatom.pdb', 'c1_t75_r1_polyAla.pdb', 'c1_t75_r1_reliable.pdb',
                 'c1_t75_r2_allatom.pdb', 'c1_t75_r2_polyAla.pdb', 'c1_t75_r2_reliable.pdb', 'c1_t75_r3_allatom.pdb',
                 'c1_t75_r3_polyAla.pdb', 'c1_t75_r3_reliable.pdb', 'c1_t80_r1_allatom.pdb', 'c1_t80_r1_polyAla.pdb',
                 'c1_t80_r1_reliable.pdb', 'c1_t80_r2_allatom.pdb', 'c1_t80_r2_polyAla.pdb', 'c1_t80_r2_reliable.pdb',
                 'c1_t80_r3_allatom.pdb', 'c1_t80_r3_polyAla.pdb', 'c1_t80_r3_reliable.pdb', 'c1_t85_r1_allatom.pdb',
                 'c1_t85_r1_polyAla.pdb', 'c1_t85_r1_reliable.pdb', 'c1_t85_r2_allatom.pdb', 'c1_t85_r2_polyAla.pdb',
                 'c1_t85_r2_reliable.pdb', 'c1_t85_r3_allatom.pdb', 'c1_t85_r3_polyAla.pdb', 'c1_t85_r3_reliable.pdb',
                 'c1_t90_r1_allatom.pdb', 'c1_t90_r1_polyAla.pdb', 'c1_t90_r1_reliable.pdb', 'c1_t90_r2_allatom.pdb',
                 'c1_t90_r2_polyAla.pdb', 'c1_t90_r2_reliable.pdb', 'c1_t90_r3_allatom.pdb', 'c1_t90_r3_polyAla.pdb',
                 'c1_t90_r3_reliable.pdb', 'c1_t95_r2_allatom.pdb', 'c1_t95_r2_polyAla.pdb', 'c1_t95_r2_reliable.pdb',
                 'c1_t95_r3_allatom.pdb', 'c1_t95_r3_polyAla.pdb', 'c1_t95_r3_reliable.pdb']

        self.assertEqual(sorted([os.path.basename(m) for m in ensembles]), eref)
        d = ensembler.ensembles_data[5]

        self.assertEqual(d['percent_truncation'], percent_truncation)
        self.assertEqual(d['truncation_method'], truncation_method)
        self.assertEqual(d['cluster_method'], cluster_method)
        self.assertEqual(d['num_clusters'], num_clusters)
        self.assertTrue(abs(d['truncation_variance'] - 13.035172) < 0001)
        self.assertEqual(d['ensemble_num_atoms'], 984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']), '4_S_00000002.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testEnsemblingThresh(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh6")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")
        
        num_clusters = 1
        cluster_method = 'spicker'
        percent_truncation = 5
        truncation_method = "thresh"
        ensembles = ensembler.generate_ensembles(models,
                                               cluster_method=cluster_method,
                                               cluster_exe=self.spicker_exe,
                                               num_clusters=num_clusters,
                                               percent_truncation=percent_truncation,
                                               truncation_method=truncation_method,
                                               work_dir=work_dir)
        
        self.assertEqual(len(ensembles), 162, len(ensembles))
        d = ensembler.ensembles_data[5]
        
        self.assertTrue(abs(d['truncation_variance'] - 27.389253) < 0001)
        self.assertEqual(d['percent_truncation'], percent_truncation)
        self.assertEqual(d['truncation_method'], truncation_method)
        self.assertEqual(d['cluster_method'], cluster_method)
        self.assertEqual(d['num_clusters'], num_clusters)
        self.assertEqual(d['subcluster_radius_threshold'], 3)
        self.assertEqual(d['side_chain_treatment'], ALLATOM)
        self.assertEqual(d['ensemble_num_atoms'], 984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']), '5_S_00000005.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testSubcluster1(self):
        """Many more then ensemble_max_models"""
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "test_subcluster")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        ensemble_max_models = 30
        
        mdir = os.path.join(self.testfiles_dir, "1p9g_models")
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'num_residues' : 5,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models_fixed_radii(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=ensemble_max_models)
        
        # Bug with theseus means cluster 1 fails
        cluster2 = [ d for d in data if d['subcluster_radius_threshold'] == 2 ][0]
        cluster3 = [ d for d in data if d['subcluster_radius_threshold'] == 3 ][0]
        
        self.assertEqual(cluster2['subcluster_num_models'], ensemble_max_models)
        self.assertEqual(cluster3['subcluster_num_models'], ensemble_max_models)
        shutil.rmtree(work_dir)
        return
    
    def testSubcluster2(self):
        """Just more then ensemble_max_models"""
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "test_subcluster")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        ensemble_max_models = 30
        
        mdir = os.path.join(self.testfiles_dir, "1p9g_models")
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")[:ensemble_max_models + 2]

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'num_residues' : 5,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models_fixed_radii(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=ensemble_max_models)
        
        # Bug with theseus means cluster 1 fails
        cluster2 = [ d for d in data if d['subcluster_radius_threshold'] == 2 ][0]
        cluster3 = [ d for d in data if d['subcluster_radius_threshold'] == 3 ][0]
        
        self.assertEqual(cluster2['subcluster_num_models'], ensemble_max_models)
        self.assertEqual(cluster3['subcluster_num_models'], ensemble_max_models)
        shutil.rmtree(work_dir)
        return

    def testSubclusterNew1(self):
        """Divergent models"""
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh7")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        ensembler.subcluster_method = "FLOATING_RADII"
        
        mdir = os.path.join(self.testfiles_dir, "2qsk_models")
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)
        
        self.assertEqual(data[0]['subcluster_num_models'], 2)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold'] - 5.82) < 0.0001)
        self.assertEqual(data[1]['subcluster_num_models'], 3)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold'] - 6.82) < 0.0001)
        self.assertEqual(data[2]['subcluster_num_models'], 4)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold'] - 6.92) < 0.0001)
        shutil.rmtree(work_dir)
        return
    
    def testSubclusterNew2(self):
        """Similar models"""
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh8")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        ensembler.subcluster_method = "FLOATING_RADII"
        
        mdir = os.path.join(self.testfiles_dir, "1p9g_models")
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)

        self.assertEqual(data[0]['subcluster_num_models'], 30,)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold'] - 1) < 0.0001)
        self.assertEqual(data[1]['subcluster_num_models'], 30)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold'] - 2) < 0.0001)
        self.assertEqual(data[2]['subcluster_num_models'], 30)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold'] - 3) < 0.0001)
        shutil.rmtree(work_dir)
        
        return
    
    def testSubclusterNew3(self):
        """standard models"""
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh9")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        ensembler.subcluster_method = "FLOATING_RADII"
        
        mdir = os.path.join(self.testfiles_dir, "models")
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)
        self.assertEqual(data[0]['subcluster_num_models'], 2)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold'] - 1.9) < 0.0001, "GOT {0}".format(data[0]['subcluster_radius_threshold']))
        self.assertEqual(data[1]['subcluster_num_models'], 3)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold'] - 2) < 0.0001, "GOT {0}".format(data[1]['subcluster_radius_threshold']))
        self.assertEqual(data[2]['subcluster_num_models'], 8)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold'] - 3) < 0.0001, "GOT {0}".format(data[2]['subcluster_radius_threshold']))
        shutil.rmtree(work_dir)
        return
    
    
    def XtestSubclusterNew(self):
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh9")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        mdir = "1P9G"
        truncated_models = glob.glob(mdir + os.sep + "*.pdb")

        subcluster_exe = self.maxcluster_exe
        clusterer = subcluster.MaxClusterer(subcluster_exe)
        clusterer.generate_distance_matrix(truncated_models)

        max_models = 30
        radius = 1
        direction = "down"
        increment = 0.1
        models, new_radius = ensembler._subcluster_nmodels(max_models,
                                                           radius,
                                                           clusterer,
                                                           direction,
                                                           increment)
        self.assertEqual(len(models), 36)
        self.assertTrue(abs(new_radius - 0.005) < 0.0001)
        return
    
    def testHomologs(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()
        ensembler.theseus_exe = self.theseus_exe
        
        work_dir = os.path.join(self.tests_dir, "homologs_test")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)

        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        alignment_file = os.path.join(self.ample_dir, 'examples', 'homologs', 'testthree.afasta')
        ensembles = ensembler.generate_ensembles_homologs(models, alignment_file=alignment_file, work_dir=work_dir)
        self.assertEqual(len(ensembles), 57)
        shutil.rmtree(work_dir)
        return
  
    def testGesamt(self):
        gesamt_exe = "/opt/ccp4-devtools/install/bin/gesamt"
        if not ample_util.is_exe(gesamt_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "gesamt_test")
        alignment_file = align_gesamt(models, gesamt_exe=gesamt_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return
    
    def testMustang(self):
        mustang_exe = "/opt/MUSTANG_v3.2.2/bin/mustang-3.2.1"
        if not ample_util.is_exe(mustang_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "mustang_test")
        alignment_file = align_mustang(models, mustang_exe=mustang_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return
    
    def testCoreFromAlignment(self):
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        work_dir = os.path.join(self.tests_dir, "homologs_core")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)

        models = glob.glob(os.path.join(self.ample_dir, "examples", "homologs", "*.pdb"))
        alignment_file = os.path.join(self.ample_dir, "examples", "homologs", "testthree.afasta")
        
        core_models = model_core_from_alignment(models, alignment_file, work_dir=work_dir)
        
        got = {}
        for m in core_models:
            name = os.path.splitext(os.path.basename(m))[0]
            resseqd = pdb_edit.resseq(m)
            resseq = resseqd[resseqd.keys()[0]]
            got[name] = resseq

        ref = { '1ujb_core':
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                56, 57, 58, 59, 60, 61, 62, 63, 64, 67, 69, 70, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85,
                86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 108, 109,
                110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130,
                131, 132, 133, 134, 135, 136, 137, 138, 139, 141, 142, 143, 144, 145, 146, 147, 149, 150, 151],
               '3c7tA_core' :
                [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154,
                 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 172, 174, 176, 180, 181, 182,
                 183, 184, 185, 186, 187, 188, 189, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 244, 245, 246, 247,
                 248, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
                 276, 277, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308,
                 309, 310, 311, 314, 315, 316],
               '2a6pA_core' :
                [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 122, 123, 124, 125, 126, 127,
                 128, 129, 130, 131, 132, 133, 134, 135, 136, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 167, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
                 179, 180, 181, 184, 185, 186, 187, 188, 189, 190, 191, 193, 194, 195] }
        
        self.assertEqual(got, ref)
        shutil.rmtree(work_dir)
        return
    
    
    def testCoreFromTheseus(self):
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        work_dir = os.path.join(self.tests_dir, "theseus_core")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        os.chdir(work_dir)

        pdbs = ['1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, "examples", "homologs", p) for p in pdbs ] 
        alignment_file = os.path.join(self.testfiles_dir, "1ujb_2a6pA_3c7tA.afasta")
        variance_file = os.path.join(self.testfiles_dir, "1ujb_2a6pA_3c7tA.variances")
        
        rt = theseus.Theseus(theseus_exe=self.theseus_exe)
        var_by_res = rt.parse_variances(variance_file)
        
        core_models = model_core_from_theseus(models, alignment_file, var_by_res, work_dir=work_dir)
        
        got = {}
        for m in core_models:
            name = os.path.splitext(os.path.basename(m))[0]
            resseqd = pdb_edit.resseq(m)
            resseq = resseqd[resseqd.keys()[0]]
            got[name] = resseq

        ref = { '1ujb_core':
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                56, 57, 58, 59, 60, 61, 62, 63, 64, 67, 69, 70, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85,
                86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 108, 109,
                110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130,
                131, 132, 133, 134, 135, 136, 137, 138, 139, 141, 142, 143, 144, 145, 146, 147, 149, 150, 151],
               '3c7tA_core' :
                [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154,
                 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 172, 174, 176, 180, 181, 182,
                 183, 184, 185, 186, 187, 188, 189, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 244, 245, 246, 247,
                 248, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
                 276, 277, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308,
                 309, 310, 311, 314, 315, 316],
               '2a6pA_core' :
                [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 122, 123, 124, 125, 126, 127,
                 128, 129, 130, 131, 132, 133, 134, 135, 136, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 167, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
                 179, 180, 181, 184, 185, 186, 187, 188, 189, 190, 191, 193, 194, 195] }
        
        self.assertEqual(got, ref)
        shutil.rmtree(work_dir)
        return

#
# Run unit tests
if __name__ == "__main__":

    if True:
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

    unittest.main(verbosity=2)

