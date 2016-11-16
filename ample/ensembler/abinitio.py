"""Ensembler module for ab initio decoys"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "17 Feb 2016"
__version__ = "1.0"

import logging
import os
import shutil

import _ensembler
import cluster_util
import subcluster
import subcluster_util
import truncation_util
from constants import SIDE_CHAIN_TREATMENTS, SUBCLUSTER_RADIUS_THRESHOLDS

from ample.util import fast_protein_cluster
from ample.util import scwrl_util
from ample.util import spicker

logger = logging.getLogger(__name__)


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
        self.subcluster_radius_thresholds = SUBCLUSTER_RADIUS_THRESHOLDS
        
        # we save the truncator so that we can query it for data later
        self.truncator = None
        
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
        logger.info('Generating {0} clusters using method: {1}'.format(num_clusters, cluster_method))

        if cluster_method != 'import' and not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
        
        # Get the cluster_method_type and cluster_score_type from the cluster_method
        cluster_method_type, cluster_score_type, cluster_exe = self.parse_cluster_method(cluster_method)
        
        # Set directory
        if cluster_method_type not in ['skip']:
            pass
        elif cluster_method_type in ['import', 'random']:
            if not os.path.isdir(cluster_dir): 
                raise RuntimeError('Cannot find cluster directory: {0}'.format(cluster_dir))
        else:
            cluster_dir = os.path.join(self.work_dir, 'clustering')
        
        if cluster_method_type == 'fast_protein_cluster':
            SCORE_TYPE = 'rmsd'
            CLUSTER_METHOD = 'kmeans'
            logger.info('Running fast_protein_cluster with: score_type: {0} cluster_method: {1}'.format(SCORE_TYPE, CLUSTER_METHOD))
            clusters = fast_protein_cluster.FPC().fpc.cluster(cluster_method=CLUSTER_METHOD,
                                                              fpc_exe=cluster_exe,
                                                              max_cluster_size=max_cluster_size,
                                                              models=models,
                                                              num_clusters=num_clusters,
                                                              nproc=self.nproc,
                                                              score_type=SCORE_TYPE,
                                                              work_dir=cluster_dir)
        elif cluster_method_type == 'import':
            clusters = cluster_util.import_cluster(cluster_dir)
        elif cluster_method_type == 'random':
            clusters = cluster_util.random_cluster(cluster_method_type,
                                                   max_cluster_size,
                                                   models,
                                                   num_clusters)
        elif cluster_method_type == 'spicker':
            logger.info('* Running SPICKER to cluster models *')
            spickerer = spicker.Spickerer(spicker_exe=cluster_exe)
            clusters = spickerer.cluster(models,
                                         num_clusters=num_clusters,
                                         max_cluster_size=max_cluster_size,
                                         score_type=cluster_score_type,
                                         run_dir=cluster_dir,
                                         score_matrix=None,
                                         nproc=self.nproc)
            logger.debug(spickerer.results_summary())
        else:
            msg = 'Unrecognised clustering method: {0}'.format(cluster_method_type)
            raise RuntimeError(msg)
                
        return clusters
    
    def ensemble_from_subcluster(self, cluster_files, radius, truncation, cluster_score=None):
        subcluster_dir = os.path.join(truncation.directory, 'subcluster_{0}'.format(radius))
        os.mkdir(subcluster_dir)
        os.chdir(subcluster_dir)

        cluster_num = truncation.cluster.index
        truncation_level = truncation.level
        basename = 'c{0}_t{1}_r{2}'.format(cluster_num, truncation_level, radius)
        
        # List of files for reference
        with open(os.path.join(subcluster_dir, "{0}.list".format(basename)), 'w') as f:
            for m in cluster_files: f.write(m + "\n")
            f.write("\n")
        
        cluster_file = self.superpose_models(cluster_files, work_dir=subcluster_dir)
        if not cluster_file:
            msg = "Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename,
                                                                                                             subcluster_dir)
            logger.critical(msg)
            return None
         
        ensemble_pdb = os.path.join(subcluster_dir, basename + '.pdb')
        shutil.move(cluster_file, ensemble_pdb)
        
        ensemble = _ensembler.Ensemble()
        
        # First add the data from the cluster
        ensemble.cluster_method = truncation.cluster.cluster_method
        ensemble.num_clusters = truncation.cluster.num_clusters
        ensemble.cluster_num = truncation.cluster.index
        ensemble.cluster_centroid = truncation.cluster.centroid
        ensemble.cluster_num_models = len(truncation.cluster)
        
        # Then the truncation information
        ensemble.num_residues = truncation.num_residues
        ensemble.truncation_dir = truncation.directory
        ensemble.truncation_level = truncation.level
        ensemble.truncation_method = truncation.method
        ensemble.truncation_percent = truncation.percent
        ensemble.truncation_residues = truncation.residues
        ensemble.truncation_variance = truncation.variances

        # Now the subcluster info
        # The data we've collected is the same for all pdbs in this level so just keep using the first
        ensemble.subcluster_num_models = len(cluster_files)
        # Get the centroid model name from the list of files given to theseus - we can't parse
        # the pdb file as theseus truncates the filename
        ensemble.subcluster_centroid_model =  os.path.abspath(cluster_files[0])
        ensemble.subcluster_radius_threshold = radius
        ensemble.subcluster_score = cluster_score
        
        ensemble.pdb = ensemble_pdb
        
        return ensemble

    def generate_ensembles(self,
                           models,
                           cluster_dir=None,
                           cluster_method=None,
                           ensembles_directory=None,
                           ensemble_max_models=None,
                           num_clusters=None,
                           percent_truncation=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS,
                           subcluster_radius_thresholds=SUBCLUSTER_RADIUS_THRESHOLDS,
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
        
        logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError, "Problem reading models given to Ensembler: {0}".format(models) 
        
        self.ensembles = []
        for cluster in self.cluster_models(models=models,
                                           cluster_method=cluster_method,
                                           num_clusters=num_clusters,
                                           cluster_dir=cluster_dir):
            if len(cluster) < 2:
                logger.info("Cannot truncate cluster {0} as < 2 models!".format(cluster.index))
                continue
            logger.info('Processing cluster: {0}'.format(cluster.index))
            
            truncate_dir = os.path.join(self.work_dir, "cluster_{0}".format(cluster.index))
            if not os.path.isdir(truncate_dir): os.mkdir(truncate_dir)
            os.chdir(truncate_dir)
            
            # Add sidechains using SCWRL here so we only add them to the models we actually use
            if use_scwrl:
                cluster.models = self.scwrl_models(cluster.models, truncate_dir, self.scwrl_exe)
                 
            self.truncator = truncation_util.Truncator(work_dir=truncate_dir)
            self.truncator.theseus_exe = self.theseus_exe
            for truncation in self.truncator.truncate_models(models=cluster.models,
                                                             truncation_method=truncation_method,
                                                             percent_truncation=percent_truncation,
                                                             truncation_pruning=truncation_pruning):
                # Add cluster information
                truncation.cluster = cluster
                for ensemble in self.subcluster_models(truncation,
                                                       subcluster_program=subcluster_program,
                                                       ensemble_max_models=self.ensemble_max_models,
                                                       radius_thresholds=subcluster_radius_thresholds):
                    # Now add the side chains
                    for ensemble in self.edit_side_chains(ensemble, side_chain_treatments):
                        self.ensembles.append(ensemble)
        return self.ensembles

    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary."""
        kwargs = {
                  'cluster_dir' : amoptd['cluster_dir'],
                  'cluster_method' : amoptd['cluster_method'],
                  'num_clusters' : amoptd['num_clusters'],
                  'percent_truncation' : amoptd['percent'],
                  'side_chain_treatments' : amoptd['side_chain_treatments'],
                  'subcluster_program' : amoptd['subcluster_program'],
                  'subcluster_radius_thresholds' : amoptd['subcluster_radius_thresholds'],
                  'truncation_method' : amoptd['truncation_method'],
                  'truncation_pruning' : amoptd['truncation_pruning'],
                  'use_scwrl' : amoptd['use_scwrl']
                  }

        # strip out any that are None
        kwargs = { k : v for k, v in kwargs.iteritems() if v is not None }
        
        ensembles = self.generate_ensembles(models, **kwargs)
        
        # We need to save these data to amopt as it's impossible to reconstruct otherwise
        amoptd['truncation_levels'] = self.truncator.truncation_levels
        amoptd['truncation_variances'] = self.truncator.truncation_variances
        amoptd['truncation_nresidues'] =  self.truncator.truncation_nresidues
        
        return ensembles

    def parse_cluster_method(self, cluster_method):
        """Return the cluster_method_type, cluster_score_type, cluster_exe from a generic cluster_method"""
        cluster_score_type = 'rmsd'
        if cluster_method == 'fast_protein_cluster':
            cluster_method_type = 'fast_protein_cluster'
            cluster_exe = self.fast_protein_cluster_exe
        elif cluster_method in ['spicker', 'spicker_qscore', 'spicker_tm']:
            cluster_method_type = 'spicker'
            cluster_exe = self.spicker_exe
            if cluster_method == 'spicker_qscore':
                cluster_score_type = 'read_matrix'
            elif cluster_method == 'spicker_tm':
                cluster_score_type = 'tm'
        elif cluster_method in ['import', 'random', 'skip']:
            cluster_method_type = cluster_method
            cluster_exe = None
        else:
            msg = "Unrecognised cluster_method: {0}".format(cluster_method)
            raise RuntimeError(msg)
        logger.debug('cluster_method_type: {0} cluster_score_type: {1} cluster_exe {2}'.format(cluster_method_type, cluster_score_type, cluster_exe))
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
                          truncation,
                          subcluster_program=None,
                          radius_thresholds=None,
                          ensemble_max_models=None):
        if self.subcluster_method == "ORIGINAL":
            return self.subcluster_models_fixed_radii(truncation,
                                                      subcluster_program=subcluster_program,
                                                      ensemble_max_models=ensemble_max_models,
                                                      radius_thresholds=radius_thresholds)
        elif self.subcluster_method == "FLOATING_RADII":
            return self.subcluster_models_floating_radii(truncation,
                                                         subcluster_program=subcluster_program,
                                                         ensemble_max_models=ensemble_max_models)
        else:
            assert False
        
    def subcluster_models_fixed_radii(self,
                                      truncation,
                                      subcluster_program=None,
                                      ensemble_max_models=None,
                                      radius_thresholds=None):
        
        # Theseus only works with > 3 residues
        if truncation.num_residues <= 2: return []
        
        if not radius_thresholds: radius_thresholds = self.subcluster_radius_thresholds
        
        # Make sure everyting happens in the truncation directory
        owd = os.getcwd()
        os.chdir(truncation.directory)
            
        # Generate the distance matrix
        clusterer = self.subclusterer_factory(subcluster_program)
        clusterer.generate_distance_matrix(truncation.models)
        # clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging

        # Loop through the radius thresholds
        ensembles = []
        previous_clusters = []
        for radius in radius_thresholds:
            logger.debug("subclustering models under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            if not cluster_files:
                logger.debug("Skipping radius {0} as no files clustered in directory {1}".format(radius, truncation.directory))
                continue
                
            logger.debug("Clustered {0} files".format(len(cluster_files)))
            cluster_files = subcluster_util.slice_subcluster(cluster_files,
                                                             previous_clusters,
                                                             ensemble_max_models,
                                                             radius,
                                                             radius_thresholds)
            if not cluster_files:
                logger.debug('Could not create different cluster for radius {0} in directory: {1}'.format(radius, truncation.directory))
                continue
            
            # Remember this cluster so we don't create duplicate clusters
            previous_clusters.append(cluster_files)
            ensemble = self.ensemble_from_subcluster(cluster_files, radius, truncation, cluster_score=clusterer.cluster_score)
            if ensemble is None: continue
            ensembles.append(ensemble)
        
        # back to where we started
        os.chdir(owd)
        
        return ensembles
    
    def subcluster_models_floating_radii(self,
                                         truncation,
                                         subcluster_program=None,
                                         ensemble_max_models=None):
        logger.info("subclustering with floating radii")

        clusterer = self.subclusterer_factory(subcluster_program)
        clusterer.generate_distance_matrix(truncation.models)
        # clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging
        
        subclusters = []
        clusters = []
        radii = []
        len_truncated_models = len(truncation.models)
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
                    logger.debug('Could not cluster files under radius: {0} - could not find different models'.format(radius))
                    break
            
            # Need to check in case we couldn't cluster under this radius
            if cluster_size == 0 or radius in radii:
                logger.debug('Could not cluster files under radius: {0} - got {1} files'.format(radius, len(cluster_files)))
                break
            logger.debug('Subclustering {0} files under radius {1}'.format(cluster_size, radius))
            cluster_ensemble = self.ensemble_from_subcluster(list(cluster_files), radius, truncation, cluster_score=clusterer.cluster_score)
            if cluster_ensemble is None:
                logger.debug('Could not cluster files under radius: {0}'.format(radius))
                break
            
            subclusters.append(cluster_ensemble)
            clusters.append(tuple(cluster_files))  # append as tuple so it is hashable
            radii.append(radius)
            if cluster_size == len_truncated_models: break
            
        return subclusters

