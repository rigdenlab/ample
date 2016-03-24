"""
02.03.2016

@author: hlfsimko
"""

import copy
import logging

_logger = logging.getLogger(__name__)

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
    
    MIN_CHUNK = 3  # We need at least 3 residues for theseus to work
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

def calculate_residues_thresh(var_by_res, percent_interval):
    """Txxx
    """

    # calculate the thresholds
    truncation_variances = generate_thresholds(var_by_res, percent_interval)

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

def generate_thresholds(var_by_res, percent_interval):
    """
    This is the original method developed by Jaclyn and used in all work until November 2014 (including the coiled-coil paper)
    
    Calculate the residue variance thresholds that will keep self.percent_interval residues for each truncation level
    """
    #--------------------------------
    # choose threshold type
    #-------------------------------
    FIXED_INTERVALS = False
    if FIXED_INTERVALS:
        thresholds = [ 1, 1.5, 2 , 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8 ]
        _logger.debug("Got {0} thresholds: {1}".format(len(thresholds), thresholds))
        return

    # List of variances ordered by residue index
    var_list = [ x.variance for x in var_by_res]
    length = len(var_list)
    if length == 0:
        msg = "Error generating thresholds, got len: {0}".format(length)
        _logger.critical(msg)
        raise RuntimeError, msg

    # How many residues should fit in each bin
    # NB - Should round up not down with int!
    chunk_size = int((float(length) / 100) * float(percent_interval))
    if chunk_size < 1:
        msg = "Error generating thresholds, got < 1 AA in chunk_size"
        _logger.critical(msg)
        raise RuntimeError, msg

    # # try to find intervals for truncation
    truncation_thresholds = _generate_thresholds(var_list, chunk_size)
    
    # Jens' new untested method
    # truncation_thresholds=self._generate_thresholds2(var_list, chunk_size)
    
    _logger.debug("Got {0} thresholds: {1}".format(len(truncation_thresholds), truncation_thresholds))

    return truncation_thresholds

def _generate_thresholds(values, chunk_size):
    """Jaclyn's threshold method
    """
    try_list = copy.deepcopy(values)
    try_list.sort()
    # print "try_list ",try_list

    # print list(chunks(try_list, int(chunk_size)))
    # For chunking list
    def chunks(a_list, chunk_size):
        for i in xrange(0, len(a_list), chunk_size):
            yield a_list[i:i + chunk_size]

    thresholds = []
    for x in list(chunks(try_list, chunk_size)):
        # print x, x[-1]
        # For some cases, multiple residues share the same variance so we don't create a separate thereshold
        if x[-1] not in thresholds:
            thresholds.append(x[-1])
            
    return thresholds

def _generate_thresholds2(values, chunk_size):
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
    """split a sequence of length into chunks each separated by percent_interval each being at least min_chunk size"""
    
    if length <= min_chunk: return [length - 1]

    # How many residues should fit in each bin
    chunk_size = int(round(float(length) * float(percent_interval) / 100.0))
    if chunk_size <= 0: return [length - 1]
    idxs = [length - 1]
    while True:
        start = idxs[-1] - chunk_size
        if start <= 0: break
        remainder = start + 1
        if remainder >= min_chunk:
            idxs.append(start)
        else:
            break
    return idxs
