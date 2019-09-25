"""Cluster utility module"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "23 Feb 2016"
__version__ = "1.0"

import logging
import random

from ample.ensembler._ensembler import Cluster

logger = logging.getLogger(__name__)


def import_cluster(cluster_models):
    """Import a cluster

    Parameters
    ----------
    cluster_models : list, tuple
       A list of models in a single cluster

    Returns
    -------
    list
       A list containing the cluster information

    """
    cluster = Cluster()
    cluster.method = "import"
    cluster.num_clusters = 1
    cluster.index = 1
    cluster.models = cluster_models
    return [cluster]


def random_cluster(cluster_method, max_cluster_size, models, num_clusters):
    """Cluster decoys using madness

    Parameters
    ----------
    cluster_method : str
       The method to be used to cluster the decoys
    max_cluster_size : int
       The maximum number of decoys per cluster
    models : list
       A list containing structure decoys
    num_clusters : int
       The number of clusters to produce

    Returns
    -------
    list
       A list containing the clusters

    Raises
    ------
    RuntimeError
       Cannot ramdonly cluster so few decoys
    RuntimeError
       Cannot find random clusters

    """
    if len(models) <= max_cluster_size + 50:  # completely arbitary number
        raise RuntimeError('Cannot randomly cluster so few models')

    i = 0
    clusters = []
    while len(clusters) < num_clusters:
        if i > num_clusters * 3:
            raise RuntimeError('Cannot find random clusters')
        cmodels = set(random.sample(models, max_cluster_size))
        if cmodels in clusters:
            logger.debug('Found duplicate cluster')
            continue

        # Data on the models
        cluster = Cluster()
        cluster.method = cluster_method
        cluster.num_clusters = num_clusters
        cluster.index = i + 1
        cluster.models = list(cmodels)  # convert list back to set

        clusters.append(cluster)
        i += 1

    return clusters
