"""Ensembler core module"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Nov 2016"
__version__ = "1.0"

import copy
import logging
import os
import re
import shutil

from constants import ENSEMBLE_MAX_MODELS, ALLATOM, POLYALA, RELIABLE, UNMODIFIED
from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import sequence_util
from ample.util import theseus

logger = logging.getLogger(__name__)


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


class Cluster(object):
    """Class to hold all data related to a cluster"""
    def __init__(self):
        self.cluster_method = None
        self.score_type = None # The scoring method used for clustering
        self.index = None # index of the cluster in the list of clusters (counting from 1)
        self.models = []
        self.num_clusters = None # How many clusters there are in list of clusters
        self.r_cen = []  # ordered list of the distance from the cluster centroid for each pdb (for Spicker)
        
        # This may be set if it's not the first model
        self._centroid = None
    
    @property
    def centroid(self):
        return self.models[0] if self._centroid is None else self._centroid
    
    @centroid.setter
    def centroid(self, model):
        self._centroid = model
    
    @property
    def size(self):
        return self.__len__()
    
    def __len__(self):
        """Return the number of models in this cluster."""
        return 0 if self.models is None else len(self.models)

    def __str__(self):
        """Return a string representation of this object."""
        _str = super(Cluster, self).__str__() + "\n"
        # Iterate through all attributes in order
        for k in sorted(self.__dict__.keys()):
            _str += "{0} : {1}\n".format(k, self.__dict__[k])
        return _str
    
class Ensemble(object):
    """Class to hold data relating to an ensemble of one or more molecular models"""
    
    def __init__(self, pdb=None):
        
        # ensemble info
        self.name = None
        self.pdb = None
        self.side_chain_treatment = None
        self.ensemble_num_atoms = None
        
        # cluster info
        self.cluster_method = None
        self.cluster_score_type = None
        self.num_clusters = None
        self.cluster_num = None
        self.cluster_centroid = None
        self.cluster_num_models = None
        
        # truncation info
        self.truncation_dir = None
        self.truncation_level = None
        self.truncation_method = None
        self.truncation_percent = None
        self.truncation_residues = None
        self.truncation_score_key = None
        self.truncation_variance = None
        self.num_residues = None
    
        # subclustering info
        self.subcluster_centroid_model = None
        self.subcluster_num_models = None
        self.subcluster_radius_threshold = None
        self.subcluster_score = None
    
        if pdb: self.from_pdb(pdb)
        return

    def copy(self):
        return copy.deepcopy(self)
    
    def __str__(self):
        """Return a string representation of this object."""
        _str = super(Ensemble, self).__str__() + "\n"
        # Iterate through all attributes in order
        for k in sorted(self.__dict__.keys()):
            _str += "{0} : {1}\n".format(k, self.__dict__[k])
        return _str


class Ensembler(object):
    """Class to generate ensembles from ab inito models (all models must have the same sequence).
    
    Attributes
    ----------
    ensembles_directory : str
        Path to the directory where the final ensembles will be stored. Should
        be outside of work_dir so it won't be deleted on purging.
    ensemble_max_models: int
        The maximum number of models that will be included in an ensemble
    nproc : int
        The number of processors that each multi-core program will be run on.
    work_dir : str
        The working directory where all the processing takes place and all intermediary
        files are kept. This may be deleted when AMPLE is run with the purge option.
        
    gesamt_exe : str
        Path to an executable
    fast_protein_cluster_exe : str
        Path to an executable
    lsqkab_exe : str
        Path to an executable
    maxcluster_exe : str
        Path to an executable
    mustang_exe : str
        Path to an executable
    scwrl_exe : str
        Path to an executable
    spicker_exe : str
        Path to an executable
    theseus_exe : str
        Path to an executable
        
    """
    def __init__(self,
                 ensembles_directory=None,
                 ensemble_max_models=ENSEMBLE_MAX_MODELS,
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
        """Set the variables required by all Ensemblers.
        
        Universal variables that are required by all Ensemblers are set on initialisation.
        These include things like the directory to store ensembles and the path to any executables
        that might be required by an Ensembler.
        
        Parameters
        ----------
        ensembles_directory : str
            Path to the directory where the final ensembles will be stored. Should
            be outside of work_dir so it won't be deleted on purging.
        ensemble_max_models: int
            The maximum number of models that will be included in an ensemble
        nproc : int
            The number of processors that each multi-core program will be run on.
        work_dir : str
            The working directory where all the processing takes place and all intermediary
            files are kept. This may be deleted when AMPLE is run with the purge option.
        gesamt_exe : str
            Path to an executable
        fast_protein_cluster_exe : str
            Path to an executable
        lsqkab_exe : str
            Path to an executable
        maxcluster_exe : str
            Path to an executable
        mustang_exe : str
            Path to an executable
        scwrl_exe : str
            Path to an executable
        spicker_exe : str
            Path to an executable
        theseus_exe : str
            Path to an executable
        **kwargs
            Arbitrary keyword arguments.
        """
        
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
        
        return
            
    def edit_side_chains(self,
                         raw_ensemble,
                         side_chain_treatments=None, 
                         homologs=False,
                         single_structure=False):
        """Mutate side chains for the supplied all-atom ensembles
        
        Parameters
        ----------
        raw_ensemble : :obj:`Ensemble`
            Ensemble object referencing the model file and holding associated data
        side_chain_treatments : :obj:`list`
            A list of the side chain treatments to be applied
        homologs : bool
            True if these ensembles have been generated from homologs
        single_structure : bool
            True if the ensemble was generated from a single structure

        Returns
        -------
        ensembles : :obj:`list`
            A `list` of :obj:`Ensemble` objects    
        """
        ensembles = []
        if side_chain_treatments is None: side_chain_treatments=[UNMODIFIED]
        for sct in side_chain_treatments:
            ensemble = raw_ensemble.copy()
            ensemble.side_chain_treatment = sct
            if homologs:
                ensemble.name = 'e{0}_{1}'.format(ensemble.truncation_level, sct)
            elif single_structure:
                ensemble.name = '{0}_t{1}_{2}'.format(ensemble.truncation_score_key, 
                                                      ensemble.truncation_level, 
                                                      sct)
            else:
                ensemble.name = 'c{0}_t{1}_r{2}_{3}'.format(ensemble.cluster_num,
                                                            ensemble.truncation_level,
                                                            ensemble.subcluster_radius_threshold,
                                                            sct)
            # create filename based on name and side chain treatment
            # fpath = ample_util.filename_append(raw_ensemble,astr=sct, directory=ensembles_directory)
            fpath = os.path.join(self.ensembles_directory, "{0}.pdb".format(ensemble.name))
            
            # Create the files
            if sct == ALLATOM or sct == UNMODIFIED:
                # For all atom just copy the file
                shutil.copy2(raw_ensemble.pdb, fpath)
            elif sct == RELIABLE:
                pdb_edit.reliable_sidechains(raw_ensemble.pdb, fpath)
            elif sct == POLYALA:
                pdb_edit.backbone(raw_ensemble.pdb, fpath)
            else:
                raise RuntimeError, "Unrecognised side_chain_treatment: {0}".format(sct)
            
            # Count the number of atoms in the ensemble-only required for benchmark mode
            natoms, nresidues = pdb_edit.num_atoms_and_residues(fpath, first=True)
            
            # Process ensemble data
            ensemble.pdb = fpath
            ensemble.ensemble_num_atoms = natoms
            # check
            assert ensemble.num_residues == nresidues, "Unmatching number of residues: {0} : {1}".format(ensemble.num_residues,
                                                                                                         nresidues)
            ensembles.append(ensemble)
                
        return ensembles

    def generate_ensembles(self, models, **kwargs):
        """Generate ensembles from models and supplied key word arguments.

        Needs to be implemented in each class.
        
        Parameters
        ----------
        models : :obj:`list`
            A list of the models to generate ensembles from
        **kwargs
            Arbitrary keyword arguments.
            
        Returns
        -------
        ensembles : :obj:`list`
            A `list` of :obj:`Ensemble` objects
        """
        raise NotImplementedError
    
    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary.
        
        Needs to be implemented in each class.
        
        This takes an AMPLE options dictionary and extracts the parameters
        from it needed to call `generate_ensembles`.

        Parameters
        ----------
        models : :obj:`list`
            A list of the models to generate ensembles from
        amoptd : :obj:`dict`
            An instance of the AMPLE options dictionary.
            
        Returns
        -------
        ensembles : :obj:`list`
            A `list` of :obj:`Ensemble` objects
        """
        raise NotImplementedError

    def superpose_models(self, models, basename='theseus', work_dir=None, homologs=False):
        run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.superpose_models(models, basename=basename, homologs=homologs)
        except Exception, e:
            logger.critical("Error running theseus: {0}".format(e))
            return False
        return run_theseus.superposed_models

