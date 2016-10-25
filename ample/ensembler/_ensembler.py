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
    def __init__(self,
                 ensembles_directory=None,
                 ensemble_max_models=30,
                 nproc=1,
                 work_dir=None,
                 # Executables
                 gesamt_exe=None,
                 fast_protein_cluster_exe=None,
                 lsqkab_exe=None,
                 maxcluster_exe=None,
                 mustang_exe=None,
                 scwrl_exe=None,
                 spicker_exe=None,
                 theseus_exe=None,
                 **kwargs
                 ):
        """Set all variables required by all ensemblers"""
        
        # For all
        self.nproc = nproc
        assert ensembles_directory and work_dir
        if not os.path.isdir(ensembles_directory): os.mkdir(ensembles_directory)
        self.ensembles_directory = ensembles_directory
        if not os.path.isdir(work_dir): os.mkdir(work_dir)
        self.work_dir = work_dir
        os.chdir(work_dir)
        
        # executables
        self.gesamt_exe = gesamt_exe
        self.fast_protein_cluster_exe = fast_protein_cluster_exe
        self.lsqkab_exe = lsqkab_exe
        self.maxcluster_exe = maxcluster_exe
        self.mustang_exe = mustang_exe
        self.scwrl_exe = scwrl_exe
        self.spicker_exe = spicker_exe
        self.theseus_exe = theseus_exe     
           
        # truncation
        self.percent_truncation = 5
        self.truncation_levels = None
        self.truncation_method = "percent"
        self.truncation_nresidues = None
        self.truncation_pruning = None
        self.truncation_scorefile = None
        self.truncation_variances = None
        
        # side chain
        
        # ensembles
        self.ensemble_max_models = ensemble_max_models
        self.ensembles = None
        self.ensembles_data = None
        
        return
            
    def edit_side_chains(self, raw_ensemble, raw_ensemble_data, side_chain_treatments, 
                         homologs=False, single_structure=False):
        """Mutate side chains for the supplied all-atom ensembles"""
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
            fpath = os.path.join(self.ensembles_directory, "{0}.pdb".format(ensemble_data['name']))
            
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

    def generate_ensembles(self, models, **kwargs):
        """Generate ensembles from models and supplied key word arguments.
        
        Needs to be implemented in each class.
        """
        raise NotImplementedError
    
    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary.
        
        Needs to be implemented in each class
        """
        raise NotImplementedError

    def superpose_models(self, models, basename=None, work_dir=None, homologs=False):
        run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.superpose_models(models, basename=basename, homologs=homologs)
        except Exception, e:
            _logger.critical("Error running theseus: {0}".format(e))
            return False
        return run_theseus.superposed_models

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
    
    @staticmethod
    def _convert_residue_scores(residue_scores):
        """Create named tuple to match store residue data"""
        scores = [ScoreVariances(idx=int(res)-1,    # Required to match Theseus
                                 resSeq=int(res),
                                 variance=float(sco)) \
                      for (res, sco) in residue_scores]
        return scores
    
    @staticmethod
    def _slice_models(data, start, end):
        """Allows us to slice a list"""
        return data[start:end]

