"""Subcluster utility module"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "02 Mar 2016"
__version__ = "1.0"

import logging
import random

logger = logging.getLogger(__name__)


def pick_nmodels(models, clusters, ensemble_max_models):
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


def slice_subcluster(cluster_files, previous_clusters, ensemble_max_models, radius, radius_thresholds):
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


def subcluster_nmodels(nmodels, radius, clusterer, direction, increment):

    MINRADIUS = 0.0001
    MAXRADIUS = 100
    
    subcluster_models = clusterer.cluster_by_radius(radius)
    len_models = len(subcluster_models) if subcluster_models else 0
    
    logger.debug("subcluster nmodels: {0} {1} {2} {3} {4}".format(len_models, nmodels, radius, direction, increment))
    if len_models == nmodels or radius < MINRADIUS or radius > MAXRADIUS:
        logger.debug("nmodels: {0} radius: {1}".format(len_models, radius))
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
        logger.debug("subcluster nmodels exceeded increment. Returning: nmodels: {0} radius: {1}".format(len(subcluster_models), radius))
        return subcluster_models, radius
        
    return subcluster_nmodels(nmodels, radius, clusterer, direction, increment)

