"""
23.02.2016

@author: hlfsimko
"""

import glob
import logging
import os
import random

# Custom
from ample.util import fast_protein_cluster
from ample.util import spicker

_logger = logging.getLogger(__name__)

#### WORST PLACE TO 'PARK' THIS FUNCTION BUT CANNOT THINK OF BETTER ONE RIGHT
#### NOW - SO ANYONE WHO SEES THIS PLEASE MOVE IF YOU FIND A BETTER PLACE
def _create_dict():
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

def fast_protein_cluster(cluster_exe, max_cluster_size, models, num_clusters, 
                         nproc, work_dir):
    """Cluster models using fast_protein_cluster program"""
    
    fpc = fast_protein_cluster.FPC()
    SCORE_TYPE = 'rmsd'
    CLUSTER_METHOD = 'kmeans'
    _logger.info('Running fast_protein_cluster with: score_type: {0} cluster_method: {1}'.format(SCORE_TYPE,
                                                                                                 CLUSTER_METHOD))
    fpc_rundir = os.path.join(work_dir, 'fast_protein_cluster')
    _logger.info('fast_protein_cluster running in directory: {0}'.format(fpc_rundir))
    
    clusters, clusters_data = fpc.cluster(cluster_method=CLUSTER_METHOD,
                                          fpc_exe=cluster_exe,
                                          max_cluster_size=max_cluster_size,
                                          models=models,
                                          num_clusters=num_clusters,
                                          nproc=nproc,
                                          score_type=SCORE_TYPE,
                                          work_dir=fpc_rundir)
    
    return clusters, clusters_data

def import_cluster(cluster_dir):
    """Import a cluster"""
    
    if not os.path.isdir(cluster_dir): 
        raise RuntimeError, "Import cluster cannot find directory: {0}".format(cluster_dir)
    
    _logger.info('Importing cluster models')
    
    cluster_models = glob.glob(os.path.join(cluster_dir, "*.pdb"))
    
    if not cluster_models: 
        raise RuntimeError, "Import cluster cannot find pdbs in directory: {0}".format(cluster_dir)
    
    # Data on the cluster
    cluster_data = _create_dict()
    cluster_data['cluster_num'] = 1
    cluster_data['cluster_centroid'] = cluster_models[0]
    cluster_data['cluster_num_models'] = len(cluster_models)
    cluster_data['cluster_method'] = "import"
    cluster_data['num_clusters'] = 1
    
    return [cluster_models], [cluster_data]

def spicker_cluster(models,
                    cluster_dir,
                    cluster_method,
                    cluster_exe,
                    num_clusters,
                    max_cluster_size,
                    cluster_score_type,
                    nproc,
                    score_matrix=None):
    """Cluster models using spicker - master function"""
            
    # Spicker Alternative for clustering
    _logger.info('* Running SPICKER to cluster models *')

    os.mkdir(cluster_dir)
    os.chdir(cluster_dir)

    spickerer = spicker.Spickerer(spicker_exe=cluster_exe)
    spickerer.cluster(models, score_type=cluster_score_type, run_dir=cluster_dir, score_matrix=score_matrix, nproc=nproc)
    _logger.debug(spickerer.results_summary())
    
    ns_clusters=len(spickerer.results)
    if ns_clusters == 0: 
        raise RuntimeError,"No clusters returned by SPICKER"
    if ns_clusters < num_clusters:
        _logger.critical('Requested {0} clusters but SPICKER only found {0} so using {1} clusters'.format(num_clusters, ns_clusters))
        num_clusters = ns_clusters
    
    clusters, clusters_data = [], []
    for i in range(num_clusters):
        # We truncate the list of models to max_cluster_size. This probably 
        # needs to be redone, because as the models are ordered by their similarity
        # to the cluster centroid, we automatically select the 200 most similar 
        # to the centroid. However if the cluster is large and the models similar
        # then when theseus calculates the variances, the variances will be 
        # representative of the 200, but might not show how the models vary throughout
        # the whole cluster, which could provide better information for truncating the models.
        #
        # Note - hlfsimko - 24.02.16: Maybe we can keep the full length clusters 
        #                             and slice the list after truncation?
        cluster = spickerer.results[i].pdbs[0:max_cluster_size]
        #cluster = spickerer.results[i].pdbs
        
        # Data on the models
        cluster_data = _create_dict()
        d = spickerer.results[i]
        cluster_data['cluster_num'] = i + 1
        cluster_data['cluster_centroid'] = d.cluster_centroid
        cluster_data['cluster_num_models'] = len(cluster)
        cluster_data['cluster_method'] = cluster_method
        cluster_data['num_clusters'] = num_clusters
        
        clusters.append(cluster)
        clusters_data.append(cluster_data)
    
    return clusters, clusters_data
    
def random_cluster(cluster_method, max_cluster_size, models, num_clusters):
    """Cluster models using madness"""
    
    clusters, clusters_data = [], []
    
    if len(models) <= max_cluster_size + 50: # completely arbitary number
        raise RuntimeError,"Cannot randomly cluster so few models!"
    
    i = 0
    while len(clusters) < num_clusters:
        if i > num_clusters * 3:
            raise RuntimeError,"Cannot find random clusters!"
        cluster = set(random.sample(models, max_cluster_size))
        if cluster in clusters:
            _logger.debug('Found duplicate cluster')
            continue
        
        # Data on the models
        cluster_data = _create_dict()
        cluster_data['cluster_num'] = i + 1
        #cluster_data['cluster_centroid'] = None
        cluster_data['cluster_num_models'] = len(cluster)
        cluster_data['cluster_method'] = cluster_method
        cluster_data['num_clusters'] = num_clusters
        
        clusters.append(cluster)
        clusters_data.append(cluster_data)
        
        i += 1
        
    # Convert all sets back to lists
    clusters = [ list(s) for s in clusters ]
    
    return clusters, clusters_data
