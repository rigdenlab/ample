#!/usr/bin/env ccp4-python

"""
17.02.2016

@author: jmht
"""

# System
import copy
import logging
import os
import shutil

# Custom
from ample.ensembler import _ensembler
from ample.ensembler import cluster_util
from ample.ensembler import subcluster
from ample.ensembler import subcluster_util
from ample.ensembler.constants import POLYALA, SIDE_CHAIN_TREATMENTS, THIN_CLUSTERS
from ample.util import scwrl_util

_logger = logging.getLogger(__name__)

class AbinitioEnsembler(_ensembler.Ensembler):
    """Ensemble creator using on multiple models with identical sequences most
       likely created using Rosetta or Quark ab initio modelling
    """
    
    def __init__(self, **kwargs):
        
        # Inherit all functions from Parent Ensembler
        super(AbinitioEnsembler, self).__init__(**kwargs)

        # self.subcluster_method='FLOATING_RADII'
        self.cluster_score_matrix = None
        self.subcluster_method = 'ORIGINAL'
        self.subcluster_program = "maxcluster"
        self.subclustering_method = "radius"
        self.subcluster_radius_thresholds = [1, 2, 3]
        
        return

    def cluster_models(self,
                       models=None,
                       cluster_method='spicker',
                       num_clusters=1,
                       cluster_dir=None,
                       max_cluster_size=200,
                       ):
        """Wrapper function to run clustering of models dependent on the method
        """
        
        # Cluster our protein structures
        _logger.info('Clustering models using method: {0}'.format(cluster_method))

        if cluster_method != 'import' and not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
        
        # Get the cluster_method_type and cluster_score_type from the cluster_method
        cluster_method_type, cluster_score_type, cluster_exe = self.parse_cluster_method(cluster_method)
        
        # Set directory
        if cluster_method_type not in ['import', 'random', 'skip']:
            cluster_dir = os.path.join(self.work_dir, 'clustering')
        
        
        if cluster_method_type == 'fast_protein_cluster':
            clusters, clusters_data = cluster_util.fast_protein_cluster(cluster_exe,
                                                                        max_cluster_size, 
                                                                        models,
                                                                        num_clusters, 
                                                                        cluster_dir)
        elif cluster_method_type == 'import':
            clusters, clusters_data = cluster_util.import_cluster(cluster_dir)
        elif cluster_method_type == 'random':
            clusters, clusters_data = cluster_util.random_cluster(cluster_method_type,
                                                                  max_cluster_size,
                                                                  models,
                                                                  num_clusters)
        elif cluster_method_type == 'spicker':
            clusters, clusters_data = cluster_util.spicker_cluster(models,
                                                                   cluster_dir,
                                                                   cluster_method_type,
                                                                   cluster_score_type,
                                                                   num_clusters,
                                                                   max_cluster_size,
                                                                   cluster_exe,
                                                                   self.nproc,
                                                                   score_matrix=self.cluster_score_matrix)
        else:
            msg = 'Unrecognised clustering method: {0}'.format(cluster_method_type)
            raise RuntimeError(msg)
                
        return clusters, clusters_data

    def generate_ensembles(self,
                           models,
                           cluster_dir=None,
                           cluster_method=None,
                           ensembles_directory=None,
                           num_clusters=None,
                           percent_truncation=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS,
                           subcluster_program=None,
                           truncation_method=None,
                           truncation_pruning=None,
                           use_scwrl=False):
 
        if not num_clusters:
            num_clusters = self.num_clusters
        if not percent_truncation:
            percent_truncation = self.percent_truncation
        if not truncation_method:
            truncation_method = self.truncation_method
        if not truncation_pruning:
            truncation_pruning = self.truncation_pruning
        
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError, "Problem reading models given to Ensembler: {0}".format(models) 
        
        _logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)

        self.ensembles = []
        self.ensembles_data = []
        for cluster_idx, (cluster_models, cluster_data) in enumerate(zip(*self.cluster_models(models=models,
                                                                                              cluster_method=cluster_method,
                                                                                              num_clusters=num_clusters,
                                                                                              cluster_dir=cluster_dir))):
            if len(cluster_models) < 2:
                _logger.info("Cannot truncate cluster {0} as < 2 models!".format(cluster_data['cluster_num']))
                continue
            _logger.info('Processing cluster: {0}'.format(cluster_idx+1))
            
            # New multi-cluster strategy
            radius_thresholds = self.subcluster_radius_thresholds
            side_chain_treatments = side_chain_treatments
            if THIN_CLUSTERS and cluster_idx > 0:
                radius_thresholds = [1, 3]
                side_chain_treatments = [ POLYALA ]
                
            truncate_dir=os.path.join(self.work_dir, "cluster_{0}".format(cluster_idx+1))
            if not os.path.isdir(truncate_dir): os.mkdir(truncate_dir)
            os.chdir(truncate_dir)
            
            # Add sidechains using SCWRL here so we only add them to the models we actually use
            if use_scwrl:
                cluster_models = self.scwrl_models(cluster_models, truncate_dir, self.scwrl_exe)
                
            for truncated_models, truncated_models_data, truncated_models_dir in zip(*self.truncate_models(cluster_models,
                                                                                                           cluster_data,
                                                                                                           truncation_method=truncation_method,
                                                                                                           truncation_pruning=truncation_pruning,
                                                                                                           percent_truncation=percent_truncation,
                                                                                                           work_dir=truncate_dir)):
                
                for subcluster, subcluster_data in zip(*self.subcluster_models(truncated_models,
                                                                               truncated_models_data,
                                                                               subcluster_program=subcluster_program,
                                                                               ensemble_max_models=self.ensemble_max_models,
                                                                               radius_thresholds=radius_thresholds,
                                                                               work_dir=truncated_models_dir)):
                    
                    for ensemble, ensemble_data in zip(*self.edit_side_chains(subcluster,
                                                                              subcluster_data,
                                                                              side_chain_treatments)):
                        self.ensembles.append(ensemble)
                        self.ensembles_data.append(ensemble_data)
        
        return self.ensembles

    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary."""
        kwargs = {'percent_truncation' : amoptd['percent'],
                  'side_chain_treatments' : amoptd['side_chain_treatments'],
                  'truncation_method' : amoptd['truncation_method'],
                  'cluster_dir' : amoptd['cluster_dir'],
                  'cluster_method' : amoptd['cluster_method'],
                  'num_clusters' : amoptd['num_clusters'],
                  'subcluster_program' : amoptd['subcluster_program'],
                  'truncation_pruning' : amoptd['truncation_pruning'],
                  'use_scwrl' : amoptd['use_scwrl']}
        return self.generate_ensembles(models, **kwargs)

    def parse_cluster_method(self, cluster_method):
        """Return the cluster_method_type, cluster_score_type, cluster_exe from a generic cluster_method"""
        cluster_score_type = 'rmsd'
        if cluster_method == 'fast_protein_cluster':
            cluster_method_type = 'fast_protein_cluster'
            cluster_exe = self.fast_protein_cluster_exe
        elif cluster_method in ['spicker', 'spicker_qm', 'spicker_tm']:
            cluster_method_type = 'spicker'
            cluster_exe = self.spicker_exe
            if cluster_method == 'spicker_qm':
                cluster_score_type = 'read_matrix'
            elif cluster_method == 'spicker_tm':
                cluster_score_type = 'tm'
        elif cluster_method in ['import', 'random', 'skip']:
            cluster_method_type = cluster_method
            cluster_exe = None
        else:
            msg = "Unrecognised cluster_method: {0}".format(cluster_method)
            raise RuntimeError(msg)
        _logger.debug('cluster_method_type: {0} cluster_score_type: {1} cluster_exe {2}'.format(cluster_method_type, cluster_score_type, cluster_exe))
        return cluster_method_type, cluster_score_type, cluster_exe
    
    def scwrl_models(self, models, work_dir, scwrl_exe):
        """Add side chains to the models with Scwrl"""
        
        scwrl_directory = os.path.join(work_dir, "scrwl")
        if not os.path.isdir(scwrl_directory): os.mkdir(scwrl_directory)
        
        scwrled_models = scwrl_util.Scwrl(scwrl_exe=scwrl_exe).process_models(models, 
                                                                              scwrl_directory, 
                                                                              strip_oxt=True)
        return scwrled_models

    def subclusterer_factory(self, subcluster_program):
        """Return an instantiated subclusterer based on the given program"""
        if subcluster_program == 'gesamt':
            clusterer = subcluster.GesamtClusterer(self.gesamt_exe, nproc=self.nproc)
        elif subcluster_program == 'maxcluster':
            clusterer = subcluster.MaxClusterer(self.maxcluster_exe)
        elif subcluster_program == 'lsqkab':
            clusterer = subcluster.LsqkabClusterer(self.lsqkab_exe)
        else:
            raise RuntimeError("Unrecognised subcluster_program: {0}".format(subcluster_program))
        return clusterer

    def subcluster_models(self,
                          truncated_models,
                          truncated_models_data,
                          subcluster_program=None,
                          radius_thresholds=None,
                          ensemble_max_models=None,
                          work_dir=None):

        if self.subcluster_method == "ORIGINAL":
            f = self.subcluster_models_fixed_radii
        elif self.subcluster_method == "FLOATING_RADII":
            f = self.subcluster_models_floating_radii
        else:
            assert False
            
        return f(truncated_models, truncated_models_data, subcluster_program, ensemble_max_models,
                 radius_thresholds=radius_thresholds, work_dir=work_dir)
        
    def subcluster_models_fixed_radii(self,
                                      truncated_models,
                                      truncated_models_data,
                                      subcluster_program=None,
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
            
        # Generate the distance matrix
        clusterer = self.subclusterer_factory(subcluster_program)
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
                                         ensemble_max_models=None,
                                         work_dir=None):
        _logger.info("subclustering with floating radii")

        clusterer = self.subclusterer_factory(subcluster_program)
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

