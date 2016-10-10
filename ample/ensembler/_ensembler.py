'''
Created on Apr 18, 2013

@author: jmht
'''

# system imports
import collections
import copy
import logging
import os
import re
import shutil

# our imports
from ample.ensembler.constants import ALLATOM, POLYALA, RELIABLE, UNMODIFIED
from ample.ensembler import cluster_util
from ample.ensembler import subcluster
from ample.ensembler import subcluster_util
from ample.ensembler import truncation_util
from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import sequence_util
from ample.util import theseus

_logger = logging.getLogger(__name__)

# Data structure to store residue information
ScoreVariances = collections.namedtuple("ScoreVariances", ["idx", "resSeq", "variance"])


def model_core_from_fasta(models, alignment_file, work_dir=None, case_sensitive=False):
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    
    # Read in alignment to get
    align_seq = sequence_util.Sequence(fasta=alignment_file)
    
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
    if case_sensitive:
        core = [ all([ x in pdb_edit.one2three.keys() for x in t ]) for t in zip(*align_seq.sequences) ]
    else:
        core = [ all([ x != GAP for x in t ]) for t in zip(*align_seq.sequences) ]

    if not any(core): raise RuntimeError("Cannot generate core for models: {0}".format(models))
    
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
    """
    Only residues from the first protein are listed in the theseus output, but then not even all of them
    
    We assume the output is based on the original alignment so that where each residue in the first protein 
    lines up with either another residue in one of the other proteins or a gap
    
    SO - we need to go through the theseus data and for each residue that is core find the corresponding residues 
    in the other proteins
    
    We use the resSeq numbers to match the residues across the alignment
    """
    if not os.path.isdir(work_dir): os.mkdir(work_dir)

    seqalign = sequence_util.Sequence(fasta=alignment_file)

    # We now need to add the list of pdbs, chains and resSeqs of the other models to the Sequence object
    for m in models: seqalign.add_pdb_data(m)
    
    # Sanity check that the names of the pdb files match those from the fasta header
    # Format is expected to be: '>1ujb.pdb(A)'
    names = [ h[1:].split('(')[0] for h in seqalign.headers ]
    if not seqalign.pdbs == names:
        raise RuntimeError, "headers and names of pdb files do not match!\n{0}\n{1}".format(seqalign.pdbs, names)
    
    # Get the name of the first pdb that the alignment is based on
    first = seqalign.pdbs[0]
    
    # Dictionary mapping model pdb to resSeqs that are core
    model2core = {}
    for p in seqalign.pdbs: model2core[p] = [] # initialise
    
    # Get list of core resSeqs in the first sequence
    model2core[first] = [ x.resSeq for x in var_by_res if x.core ]
    
    # Now go through the first sequence and get the resSeqs of the corresponding core for the other models
    pointer = 0 # Tracks where we are in the first sequence
    for i, resSeq in enumerate(seqalign.resseqs[0]):
        if model2core[first][pointer] == resSeq:
            # Core residue in first sequence so append the corresponding resSeqs for the other proteins
            for j, pdb in enumerate(seqalign.pdbs[1:]):
                model2core[pdb].append(seqalign.resseqs[j+1][i])
            pointer += 1
            if pointer >= len(model2core[first]): break
            
    core_models = []
    for m in models:
        name = os.path.basename(m)
        pdbout = ample_util.filename_append(m, astr='core', directory=work_dir)
        pdb_edit.select_residues(m, pdbout, tokeep=model2core[name])
        core_models.append(pdbout)
        
    return core_models

class Ensembler(object):
    """Class to generate ensembles from ab inito models (all models must have same sequence)
    
    """
    def __init__(self):
        
        # For all
        self.work_dir = None  # top directory where everything gets done
        self.theseus_exe = None
        self.gesamt_exe = None
        self.lsqkab_exe = os.path.join(os.environ['CCP4'],'bin','lsqkab' + ample_util.EXE_EXT)
        assert ample_util.is_exe(self.lsqkab_exe),"Cannot find lsqkab: {0}".format(self.lsqkab_exe)
        
        # clustering
        self.cluster_method = "spicker"  # the method for initial clustering
        self.cluster_exe = None
        self.num_clusters = 1
        
        # truncation
        self.percent_truncation = 5
        self.truncation_levels = None
        self.truncation_method = "percent"
        self.truncation_nresidues = None
        self.truncation_pruning = None
        self.truncation_scorefile = None
        self.truncation_variances = None
        
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
        
        self.score_matrix = None
        
        return
    
    def superpose_models(self, models, basename=None, work_dir=None, homologs=False):
        run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.superpose_models(models, basename=basename, homologs=homologs)
        except Exception, e:
            _logger.critical("Error running theseus: {0}".format(e))
            return False
        return run_theseus.superposed_models
            
    def cluster_models(self, cluster_dir=None, cluster_exe=None, cluster_method=None,
                       models=None, max_cluster_size=200, num_clusters=None, nproc=1):
        """Wrapper function to run clustering of models dependent on the method
        """
        ########################################################################
        ## Switch statement solution to find function to use
        ## `Key` equal to cluster_method
        ## `Value` equal to function name, function arguments
        switch = {'fast_protein_cluster' : [cluster_util.fast_protein_cluster, 
                                            list([cluster_exe, max_cluster_size, 
                                                  models, num_clusters, nproc, 
                                                  self.work_dir])],
                  
                  'import' : [cluster_util.import_cluster, 
                              list([cluster_dir])],
                  
                  'random' : [cluster_util.random_cluster, 
                              list([cluster_method, max_cluster_size, models, 
                                    num_clusters])],
                  
                  'spicker' : [cluster_util.spicker_default, 
                               list([cluster_exe, cluster_method, max_cluster_size, num_clusters,
                                     models, self.work_dir, nproc])],
                  
                  'spicker_qscore' : [cluster_util.spicker_qscore, 
                                      list([cluster_exe, cluster_method, max_cluster_size, num_clusters, 
                                            models, self.work_dir, nproc, self.gesamt_exe])],
                  
                  'spicker_tmscore' : [cluster_util.spicker_tmscore, 
                                       list([cluster_exe, cluster_method, max_cluster_size, num_clusters, 
                                             models, self.work_dir, nproc])],
        }
        
        # Get the function handler and keyword arguments
        cluster_function, args = switch.get(cluster_method, None)
        if not cluster_function:
            msg = 'Unrecognised clustering method: {0}'.format(cluster_method)
            raise RuntimeError(msg)
                
        # Cluster our protein structures
        _logger.info('Clustering models using method: {0}'.format(cluster_method))
        clusters, clusters_data = cluster_function(*args)
        
        return clusters, clusters_data
    
    def edit_side_chains(self, raw_ensemble, raw_ensemble_data, side_chain_treatments, 
                         ensembles_directory, homologs=False, single_structure=False):
        
        assert os.path.isdir(ensembles_directory), "Cannot find ensembles directory: {0}".format(ensembles_directory)
        ensembles = []
        ensembles_data = []
        if side_chain_treatments is None: side_chain_treatments=[UNMODIFIED]
        for sct in side_chain_treatments:
            ensemble_data = copy.copy(raw_ensemble_data)
            ensemble_data['side_chain_treatment'] = sct
            if homologs:
                ensemble_data['name'] = 'e{0}_{1}'.format(ensemble_data['truncation_level'], sct)
            elif single_structure:
                ensemble_data['name'] = '{0}_t{1}_{2}'.format(ensemble_data['truncation_score_key'], 
                                                              ensemble_data['truncation_level'], 
                                                              sct)
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

    def subcluster_models(self, truncated_models, truncated_models_data,
                          subcluster_program=None, subcluster_exe=None,
                          radius_thresholds=None, ensemble_max_models=None,
                          work_dir=None):
        
        if self.subcluster_method == "ORIGINAL":
            f = self.subcluster_models_fixed_radii
        elif self.subcluster_method == "FLOATING_RADII":
            f = self.subcluster_models_floating_radii
        else:
            assert False
            
        return f(truncated_models, truncated_models_data, subcluster_program, 
                 subcluster_exe, ensemble_max_models, radius_thresholds=radius_thresholds, 
                 work_dir=work_dir)
        
    def subcluster_models_fixed_radii(self,
                                      truncated_models,
                                      truncated_models_data,
                                      subcluster_program=None,
                                      subcluster_exe=None,
                                      ensemble_max_models=None,
                                      radius_thresholds=None,
                                      work_dir=None):
        
        # Theseus only works with > 3 residues
        if truncated_models_data['num_residues'] <= 2: return [], []
        
        if not radius_thresholds: radius_thresholds = self.subcluster_radius_thresholds
        ensembles = []
        ensembles_data = []
        
        # Use first model to get data on level
        cluster_num = truncated_models_data['cluster_num']
        truncation_level = truncated_models_data['truncation_level']
        
        # Make sure everyting happens in the truncation directory
        owd = os.getcwd()
        os.chdir(work_dir)
            
        # Run maxcluster to generate the distance matrix
        if subcluster_program == 'maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        elif subcluster_program == 'lsqkab':
            clusterer = subcluster.LsqkabClusterer(self.subcluster_exe)
        else:
            assert False,subcluster_program
        clusterer.generate_distance_matrix(truncated_models)
        # clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging

        # Loop through the radius thresholds
        previous_clusters = []
        for radius in radius_thresholds:
            _logger.debug("subclustering models under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            if not cluster_files:
                _logger.debug("Skipping radius {0} as no files clustered in directory {1}".format(radius, work_dir))
                continue
                
            _logger.debug("Clustered {0} files".format(len(cluster_files)))
            cluster_files = subcluster_util.slice_subcluster(cluster_files, previous_clusters, ensemble_max_models, radius, radius_thresholds)
            if not cluster_files:
                _logger.debug('Could not create different cluster for radius {0} in directory: {1}'.format(radius, work_dir))
                continue
            
            # Remember this cluster so we don't create duplicate clusters
            previous_clusters.append(cluster_files)

            # Got files so create the directories
            subcluster_dir = os.path.join(work_dir, 'subcluster_{0}'.format(radius))
            os.mkdir(subcluster_dir)
            os.chdir(subcluster_dir)
            basename = 'c{0}_t{1}_r{2}'.format(cluster_num, truncation_level, radius)
            
            # List of files for reference
            with open(os.path.join(subcluster_dir, "{0}.list".format(basename)), 'w') as f:
                for m in cluster_files: f.write(m + "\n")
                f.write("\n")
            
            cluster_file = self.superpose_models(cluster_files, work_dir=subcluster_dir)
            if not cluster_file:
                msg = "Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename, subcluster_dir)
                _logger.critical(msg)
                continue
             
            ensemble = os.path.join(subcluster_dir, basename + '.pdb')
            shutil.move(cluster_file, ensemble)

            # The data we've collected is the same for all pdbs in this level so just keep using the first  
            ensemble_data = copy.copy(truncated_models_data)
            ensemble_data['subcluster_num_models'] = len(cluster_files)
            ensemble_data['subcluster_radius_threshold'] = radius
            ensemble_data['subcluster_score'] =  clusterer.cluster_score
            ensemble_data['ensemble_pdb'] = ensemble

            # Get the centroid model name from the list of files given to theseus - we can't parse
            # the pdb file as theseus truncates the filename
            ensemble_data['subcluster_centroid_model'] = os.path.abspath(cluster_files[0])
            
            ensembles.append(ensemble)
            ensembles_data.append(ensemble_data)
        
        # back to where we started
        os.chdir(owd)
        
        return ensembles, ensembles_data
    
    def subcluster_models_floating_radii(self,
                                         truncated_models,
                                         truncated_models_data,
                                         subcluster_program=None,
                                         subcluster_exe=None,
                                         ensemble_max_models=None,
                                         work_dir=None):
        _logger.info("subclustering with floating radii")

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
                cluster_files, radius = subcluster_util.subcluster_nmodels(nmodels, radius, clusterer, direction='up', increment=1)
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
                cluster_files, radius = subcluster_util.subcluster_nmodels(nmodels, radius, clusterer, direction='up', increment=1)
                cluster_files = sorted(cluster_files)
            elif cluster_size >= ensemble_max_models or cluster_files in clusters:
                # Randomly pick ensemble_max_models
                cluster_files = subcluster_util.pick_nmodels(cluster_files, clusters, ensemble_max_models)
                if not cluster_files:
                    _logger.debug('Could not cluster files under radius: {0} - could not find different models'.format(radius))
                    break
            
            # Need to check in case we couldn't cluster under this radius
            if cluster_size == 0 or radius in radii:
                _logger.debug('Could not cluster files under radius: {0} - got {1} files'.format(radius, len(cluster_files)))
                break
            
            _logger.debug('Subclustering {0} files under radius {1}'.format(cluster_size, radius))
            try:
                cluster_ensemble, data = subcluster_util.subcluster_radius(list(cluster_files), radius, truncated_models_data)
            except RuntimeError:
                _logger.debug('Could not cluster files under radius: {0} as theseus failed'.format(radius, len(cluster_files)))
                # If theseus fails, we just move
                break
            
            subclusters.append(cluster_ensemble)
            subclusters_data.append(data)
            clusters.append(tuple(cluster_files))  # append as tuple so it is hashable
            radii.append(radius)
            if cluster_size == len_truncated_models: break
            
        return subclusters, subclusters_data

    def truncate_models(self,
                        models,
                        models_data={},
                        max_cluster_size=200,
                        truncation_method=None,
                        percent_truncation=None,
                        truncation_pruning=None,
                        residue_scores=None,
                        homologs=False,
                        alignment_file=None,
                        work_dir=None):
        
        assert (len(models) > 1 or residue_scores), "Cannot truncate as < 2 models!"
        assert truncation_method and percent_truncation, "Missing arguments: {0}".format(truncation_method)

        # Create the directories we'll be working in
        assert work_dir and os.path.isdir(work_dir), "truncate_models needs a work_dir"
        os.chdir(work_dir)
        
        # Calculate variances between pdb and align them (we currently only require the aligned models for homologs)
        if truncation_method != "scores":
            run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
            try: run_theseus.superpose_models(models, homologs=homologs, alignment_file=alignment_file)
            except RuntimeError as e:
                _logger.critical(e)
                return [],[]
        
        if homologs:
            # If using homologs, now trim down to the core. We only do this here so that we are using the aligned models from
            # theseus, which makes it easier to see what the truncation is doing.
            models = model_core_from_fasta(run_theseus.aligned_models,
                                           alignment_file=alignment_file,
                                           work_dir=os.path.join(work_dir,'core_models'))
            # Unfortunately Theseus doesn't print all residues in its output format, so we can't use the variances we calculated before and
            # need to calculate the variances of the core models 
            try: run_theseus.superpose_models(models, homologs=homologs, basename='homologs_core')
            except RuntimeError as e:
                _logger.critical(e)
                return [],[]
        
        # No THESEUS variances required if scores for each residue provided
        var_by_res = run_theseus.var_by_res() if truncation_method != "scores" \
            else self._convert_residue_scores(residue_scores)
            
        if not len(var_by_res) > 0:
            msg = "Error reading residue variances!"
            _logger.critical(msg)
            raise RuntimeError(msg)
        
        _logger.info('Using truncation method: {0}'.format(truncation_method))
        # Calculate which residues to keep under the different methods
        truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = None, None, None, None
        if truncation_method == 'percent':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = truncation_util.calculate_residues_percent(var_by_res, percent_truncation)
        elif truncation_method == 'scores':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = truncation_util.calculate_residues_percent(var_by_res, percent_truncation)
        elif truncation_method == 'thresh':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = truncation_util.calculate_residues_thresh(var_by_res, percent_truncation)
        elif truncation_method == 'focussed':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = truncation_util.calculate_residues_focussed(var_by_res)
        else:
            raise RuntimeError, "Unrecognised ensembling mode: {0}".format(truncation_method)
        
        self.truncation_levels = truncation_levels # save so we can put in results dict
        self.truncation_variances = truncation_variances # save so we can put in results dict
        self.truncation_nresidues = [len(r) for r in truncation_residues] # save so we can put in results dict
        
        # Use all models in cluster to calculate variance but then slice to max_cluster_size
        #models = self._slice_models(models, 0, max_cluster_size)
        
        truncated_models = []
        truncated_models_data = []
        truncated_models_dirs = []
        pruned_residues = None
        for tlevel, tvar, tresidues, tresidue_idxs in zip(truncation_levels, 
                                                          truncation_variances, 
                                                          truncation_residues, 
                                                          truncation_residue_idxs):
            # Prune singletone/doubletone etc. residues if required
            _logger.debug("truncation_pruning: {0}".format(truncation_pruning))
            if truncation_pruning == 'single':
                tresidue_idxs, pruned_residues=truncation_util.prune_residues(tresidue_idxs, chunk_size=1, allowed_gap=2)
                if pruned_residues: _logger.debug("prune_residues removing: {0}".format(pruned_residues))
            elif truncation_pruning is None:
                pass
            else:
                raise RuntimeError("Unrecognised truncation_pruning: {0}".format(truncation_pruning))
            
            # Skip if there are no residues
            if not tresidue_idxs:
                _logger.debug("Skipping truncation level {0} with variance {1} as no residues".format(tlevel, tvar))
                continue
            
            trunc_dir = os.path.join(work_dir, 'tlevel_{0}'.format(tlevel))
            os.mkdir(trunc_dir)
            _logger.info('Truncating at: {0} in directory {1}'.format(tlevel, trunc_dir))
            
            # list of models for this truncation level
            level_models = []
            for infile in models:
                pdbout = ample_util.filename_append(infile, str(tlevel), directory=trunc_dir)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdb_edit.select_residues(pdbin=infile, pdbout=pdbout, tokeep_idx=tresidue_idxs)
                level_models.append(pdbout)
            
            # Add the model
            truncated_models.append(level_models)
            truncated_models_dirs.append(trunc_dir)

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
            
        return truncated_models, truncated_models_data, truncated_models_dirs

    def _convert_residue_scores(self, residue_scores):
        """Create named tuple to match store residue data"""
        scores = [ScoreVariances(idx=int(res)-1,    # Required to match Theseus
                                 resSeq=int(res),
                                 variance=float(sco)) \
                      for (res, sco) in residue_scores]
        return scores
    
    def _slice_models(self, data, start, end):
        """Allows us to slice a list"""
        return data[start:end]

