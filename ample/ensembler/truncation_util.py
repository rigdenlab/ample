"""Truncation utility module"""
__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "10 Jul 2018"
__version__ = "1.0"

import collections
from enum import Enum
import logging
import os

from ample.ensembler._ensembler import model_core_from_fasta 
from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import theseus

logger = logging.getLogger(__name__)

MIN_CHUNK = 3 # Theseus needs at least 3 residues in order to work

class TRUNCATION_METHODS(Enum):
    FOCUSED = 'focussed'
    PERCENT = 'percent'
    PERCENT_FIXED = 'percent_fixed_intervals'
    SCORES = 'scores'

# Data structure to store residue information
ScoreVariances = collections.namedtuple("ScoreVariances", ["idx", "resSeq", "variance"])


def calculate_residues_focussed(var_by_res):
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
        return calculate_residues_percent(var_by_res, 5)

    # Get list of residue indices sorted by variance - from least variable to most
    var_by_res.sort(key=lambda x: x.variance, reverse=False)
     
    # Split a 40 - length interval into 10 even chunks.
    llen = 40
    lower_start = _split_sequence(llen, 10)
    
    # Split remaining interval into 10 even chunks. We need to add the start sequence as we have
    # removed llen residues
    ulen = length - llen
    upper_start = [ i + llen for i in _split_sequence(ulen, 10) ]
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


def calculate_residues_percent(var_by_res, percent_interval):
    """Calculate the list of residues to keep if we are keeping self.percent residues under
    each truncation bin. The threshold is just the threshold of the most variable residue"""
    
    length = len(var_by_res)
    start_idxs = _split_sequence(length, percent_interval, min_chunk=MIN_CHUNK)
    
    # Get list of residue indices sorted by variance - from least to most
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
        idxs = idxs_all[:start + 1]
        thresh = variances[start]  # For the threshold we take the threshold of the most variable residue
        truncation_variances.append(thresh)
        truncation_levels.append(percent)
        truncation_residues.append(sorted(residues))
        truncation_residue_idxs.append(sorted(idxs))
             
    return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs


def calculate_residues_percent_fixed_intervals(var_by_res, percent_fixed_intervals):
    """Calculate the list of residues to keep if we are keeping the percentage of residues specified
    in the list of integers percent_fixed_intervals"""
    
    if not (isinstance(percent_fixed_intervals, list) and \
           len(percent_fixed_intervals) > 0 and \
            min(percent_fixed_intervals) > 0 and \
            max(percent_fixed_intervals) <= 100):
        raise RuntimeError("Incorrect percent_fixed_intervals argument: {}".format(percent_fixed_intervals))
    
    # Get list of residue indices sorted by variance - from least to most
    var_by_res.sort(key=lambda x: x.variance, reverse=False)
    idxs_all = [x.idx for x in var_by_res] # indexes correspond to the number of the residue - 1
    resseq_all = [x.resSeq for x in var_by_res]
    variances = [x.variance for x in var_by_res]
    length = len(var_by_res)
    
    def find_closest_index(nresidues_ideal, idxs_all):
        min_dist = 99999 
        start = None
        for idx in idxs_all:
            rdist = abs(nresidues_ideal - (idx + 1))
            if  rdist <= min_dist and idx + 1 >= MIN_CHUNK:
                min_dist = rdist
                start = idx
        assert start is not None
        return start
    
    # Get list of residues to keep under the different intevals
    truncation_levels = []
    truncation_variances = []
    truncation_residues = []
    truncation_residue_idxs = []
    for interval in percent_fixed_intervals:
        nresidues_ideal = (length * int(interval)) / 100
        start = find_closest_index(nresidues_ideal, idxs_all)
        percent = int(round(float(start + 1) / float(length) * 100))
        residues = resseq_all[:start + 1]
        idxs = resseq_all[:start + 1]
        idxs = idxs_all[:start + 1]
        thresh = variances[start]  # For the threshold we take the threshold of the most variable residue
        truncation_variances.append(thresh)
        truncation_levels.append(percent)
        truncation_residues.append(sorted(residues))
        truncation_residue_idxs.append(sorted(idxs))
             
    return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs


def prune_residues(residues, chunk_size=1, allowed_gap=2):
    """Remove any residues that are < chunk_size where the gap before and after is > allowed_gap"""
    
    assert chunk_size > 0 and allowed_gap > 0, \
        "chunk_size and allowed_gap must be > 0!: {0} {1}".format(chunk_size, allowed_gap)
    
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
    
    idxLast = lenr - 1
    for i in xrange(1, idxLast+1):
        this_residue = residues[i]
        
        if i == idxLast or this_residue != last + 1:
        
            if i == idxLast and this_residue != last + 1:
                start = this_residue
                last_chunk_end = last
                last = this_residue
                postgap = allowed_gap + 1
        
            elif i == idxLast and this_residue == last + 1:
                last = this_residue
                postgap = allowed_gap + 1
                  
            elif i != idxLast and this_residue != last + 1:
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
            
        last = this_residue
    
    # Remove the chunks and return
    if len(to_remove):
        return [r for r in residues if r not in to_remove], to_remove
    else:
        return residues, None


def _split_sequence(length, percent_interval, min_chunk=3):
    """Split a sequence of length into chunks each separated by percent_interval each being at least min_chunk size"""
    if length <= min_chunk:
        return [length - 1]
    # How many residues should fit in each bin?
    chunk_size = int(round(float(length) * float(percent_interval) / 100.0))
    if chunk_size <= 0:
        return [length - 1]
    idxs = [length - 1]
    while True:
        start = idxs[-1] - chunk_size
        if start <= 0:
            break
        remainder = start + 1
        if remainder >= min_chunk:
            idxs.append(start)
        else:
            break
    return idxs


class Truncation(object):
    """Holds information relating to a single truncation of a cluster of models"""
    def __init__(self):
        self.cluster = None # The cluster object this truncation was created from
        self.directory = None
        self.level = None
        self.method = None
        self.models = None
        self.percent = None
        self.residues = None
        self.residues_idxs = None
        self.variances = None
    
    @property
    def num_residues(self):
        return 0 if self.residues is None else len(self.residues)

    def __str__(self):
        """Return a string representation of this object."""
        _str = super(Truncation, self).__str__() + "\n"
        # Iterate through all attributes in order
        for k in sorted(self.__dict__.keys()):
            _str += "{0} : {1}\n".format(k, self.__dict__[k])
        return _str
        

class Truncator(object):
    def __init__(self, work_dir):
        """Class to take one or more models and truncate them based on a supplied or generated metric"""
        self.work_dir = work_dir
        self.models = None
        self.aligned_models = None
        self.truncations = None
        self.theseus_exe = None
        
        # We keep these for bookeeping as they go in the ample dictionary
        self.truncation_levels =  None
        self.truncation_variances = None
        self.truncation_nresidues = None
    
    def calculate_truncations(self,
                              models=None,
                              truncation_method=None,
                              percent_truncation=None,
                              percent_fixed_intervals=None,
                              truncation_pruning=None,
                              residue_scores=None,
                              alignment_file=None,
                              homologs=False):
        """Returns a list of Truncation objects, one for each truncation level.
        
        This method doesn't do any truncating - it just calculates the data for each truncation level.
        """
        
        assert (len(models) > 1 or residue_scores), "Cannot truncate as < 2 models!"
        assert truncation_method and percent_truncation, "Missing arguments: {0} : {1}".format(truncation_method,
                                                                                               percent_truncation)
        assert ample_util.is_exe(self.theseus_exe),"Cannot find theseus_exe: {0}".format(self.theseus_exe)

        # Create the directories we'll be working in
        assert self.work_dir and os.path.isdir(self.work_dir), "truncate_models needs a self.work_dir"
        os.chdir(self.work_dir)
        
        self.models = models
        # Calculate variances between pdb and align them (we currently only require the aligned models for homologs)
        if truncation_method != "scores":
            run_theseus = theseus.Theseus(work_dir=self.work_dir, theseus_exe=self.theseus_exe)
            try:
                run_theseus.superpose_models(self.models, homologs=homologs, alignment_file=alignment_file)
                self.aligned_models = run_theseus.aligned_models
            except RuntimeError as e:
                logger.critical(e)
                return []
        
        if homologs:
            # If using homologs, now trim down to the core. We only do this here so that we are using the aligned models from
            # theseus, which makes it easier to see what the truncation is doing.
            models = model_core_from_fasta(self.aligned_models,
                                           alignment_file=alignment_file,
                                           work_dir=os.path.join(self.work_dir,'core_models'))
            # Unfortunately Theseus doesn't print all residues in its output format, so we can't use the variances we calculated before and
            # need to calculate the variances of the core models 
            try:
                run_theseus.superpose_models(models, homologs=homologs, basename='homologs_core')
                self.models = run_theseus.aligned_models
                self.aligned_models = run_theseus.aligned_models
            except RuntimeError as e:
                logger.critical(e)
                return []
        
        # No THESEUS variances required if scores for each residue provided
        var_by_res = run_theseus.var_by_res if truncation_method != "scores" \
            else self._convert_residue_scores(residue_scores)
        if not len(var_by_res) > 0:
            raise RuntimeError("Error reading residue variances!")
        
        logger.info('Using truncation method: %s', truncation_method)
        # Calculate which residues to keep under the different methods
        if truncation_method in [TRUNCATION_METHODS.PERCENT, TRUNCATION_METHODS.SCORES]:
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = calculate_residues_percent(var_by_res, percent_truncation)
        elif truncation_method == TRUNCATION_METHODS.PERCENT_FIXED:
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = calculate_residues_percent_fixed_intervals(var_by_res, percent_fixed_intervals)
        elif truncation_method == TRUNCATION_METHODS.FOCUSED:
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = calculate_residues_focussed(var_by_res)
        else:
            raise RuntimeError("Unrecognised ensembling mode: {0}".format(truncation_method))
        
        # Somewhat of a hack to save the data so we can put it in the amoptd
        self.truncation_levels =  truncation_levels
        self.truncation_variances = truncation_variances
        self.truncation_nresidues =  [len(r) for r in truncation_residues]
        
        truncations = []
        for tlevel, tvar, tresidues, tresidue_idxs in zip(truncation_levels, 
                                                          truncation_variances, 
                                                          truncation_residues, 
                                                          truncation_residue_idxs):
            # Prune singletone/doubletone etc. residues if required
            logger.debug("truncation_pruning: %s", truncation_pruning)
            if truncation_pruning == 'single':
                tresidue_idxs, pruned_residues = prune_residues(tresidue_idxs, chunk_size=1, allowed_gap=2)
                if pruned_residues:
                    logger.debug("prune_residues removing: %s", pruned_residues)
            elif truncation_pruning is None:
                pass
            else:
                raise RuntimeError("Unrecognised truncation_pruning: {0}".format(truncation_pruning))
            
            # Skip if there are no residues
            if not tresidue_idxs:
                logger.debug("Skipping truncation level %s with variance %s as no residues", tlevel, tvar)
                continue
            truncation = Truncation()
            truncation.method = truncation_method
            truncation.percent = percent_truncation
            truncation.level = tlevel
            truncation.variances = tvar
            truncation.residues = tresidues
            truncation.residues_idxs = tresidue_idxs
            truncations.append(truncation)
        return truncations
      
    def truncate_models(self,
                        models,
                        max_cluster_size=200,
                        truncation_method=None,
                        percent_truncation=None,
                        percent_fixed_intervals=None,
                        truncation_pruning=None,
                        residue_scores=None,
                        homologs=False,
                        alignment_file=None,
                        work_dir=None):
        """Generate a set of Truncation objects, referencing a set of truncated models generated from the supplied models"""
        truncations = self.calculate_truncations(models=models,
                                                 truncation_method=truncation_method,
                                                 percent_truncation=percent_truncation,
                                                 percent_fixed_intervals=percent_fixed_intervals,
                                                 truncation_pruning=truncation_pruning,
                                                 residue_scores=residue_scores,
                                                 alignment_file=alignment_file,
                                                 homologs=homologs)
        if truncations is None or len(truncations) < 1:
            logger.critical("Unable to truncate the ensembles - no viable truncations")
            return []
        # Loop through the Truncation objects, truncating the models based on the truncation data and adding
        # the truncated models to the Truncation.models attribute
        for truncation in truncations:
            truncation.directory = os.path.join(self.work_dir, 'tlevel_{0}'.format(truncation.level))
            os.mkdir(truncation.directory)
            logger.info('Truncating at: %s in directory %s', truncation.level, truncation.directory)
            truncation.models = []
            for infile in self.models:
                pdbout = ample_util.filename_append(infile, str(truncation.level), directory=truncation.directory)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdb_edit.select_residues(pdbin=infile, pdbout=pdbout, tokeep_idx=truncation.residues_idxs)
                truncation.models.append(pdbout)
        self.truncations = truncations
        return truncations

    @staticmethod
    def _convert_residue_scores(residue_scores):
        """Create named tuple to match store residue data"""
        scores = [ScoreVariances(idx=int(res)-1,    # Required to match Theseus
                                 resSeq=int(res),
                                 variance=float(sco)) \
                      for (res, sco) in residue_scores]
        return scores

