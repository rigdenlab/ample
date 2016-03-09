#!/usr/bin/env ccp4-python

"""
17.02.2016

@author: jmht
"""

# System
import glob
import logging
import os
import shutil
import sys
import unittest

# Custom
import ample_scwrl
import ample_util
import ensembler

# Inherit some variables defined in the ensembler module
ALLATOM = ensembler.ALLATOM
POLYALA = ensembler.POLYALA
RELIABLE = ensembler.RELIABLE
SIDE_CHAIN_TREATMENTS = ensembler.SIDE_CHAIN_TREATMENTS
THIN_CLUSTERS = ensembler.THIN_CLUSTERS
UNMODIFIED = ensembler.UNMODIFIED

_logger = logging.getLogger(__name__)


class Ensembler(ensembler.Ensembler):
    """Ensemble creator using on multiple models with identical sequences most
       likely created using Rosetta or Quark ab initio modelling
    """
    
    def __init__(self):
        # Inherit all functions from Parent Ensembler
        ensembler.Ensembler.__init__(self)
        
        self.scwrl_exe = None
        
        return
    
    def generate_ensembles(self,
                           models,
                           cluster_dir=None,
                           cluster_exe=None,
                           cluster_method=None,
                           ensembles_directory=None,
                           nproc=None,
                           num_clusters=None,
                           percent_truncation=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS,
                           truncation_method=None,
                           truncation_pruning=None,
                           use_scwrl=False,
                           work_dir=None):
        
        # Work dir set each time
        if not work_dir: raise RuntimeError, "Need to set work_dir!"
        self.work_dir = work_dir
        
        if not cluster_method:
            cluster_method = self.cluster_method
        if not cluster_exe:
            cluster_exe = self.cluster_exe
        if not num_clusters:
            num_clusters = self.num_clusters
        if not percent_truncation:
            percent_truncation = self.percent_truncation
        if not truncation_method:
            truncation_method = self.truncation_method
        if not truncation_pruning:
            truncation_pruning = self.truncation_pruning
        if not ensembles_directory:
            self.ensembles_directory = os.path.join(work_dir, "ensembles")
        else:
            self.ensembles_directory = ensembles_directory
        
        if not cluster_method is 'import' and not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
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
                                                                                              cluster_exe=cluster_exe,
                                                                                              cluster_dir=cluster_dir,
                                                                                              nproc=nproc))):
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
                                                                               subcluster_program=self.subcluster_program,
                                                                               subcluster_exe=self.subcluster_program,
                                                                               ensemble_max_models=self.ensemble_max_models,
                                                                               radius_thresholds=radius_thresholds,
                                                                               work_dir=truncated_models_dir)):
                    
                    for ensemble, ensemble_data in zip(*self.edit_side_chains(subcluster,
                                                                              subcluster_data,
                                                                              side_chain_treatments,
                                                                              self.ensembles_directory)):
                        self.ensembles.append(ensemble)
                        self.ensembles_data.append(ensemble_data)
        
        return self.ensembles
    
    def scwrl_models(self, models, work_dir, scwrl_exe):
        """Use Scwrl to add side chains to the models"""
        
        scwrl_directory = os.path.join(work_dir, "scrwl")
        if not os.path.isdir(scwrl_directory): os.mkdir(scwrl_directory)
        
        scwrled_models = ample_scwrl.Scwrl(scwrl_exe=scwrl_exe).process_models(models, 
                                                                               scwrl_directory, 
                                                                               strip_oxt=True)
        return scwrled_models
        

