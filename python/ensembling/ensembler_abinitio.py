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
        

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        
        cls.theseus_exe = ample_util.find_exe('theseus')
        cls.spicker_exe = ample_util.find_exe('spicker')
        cls.maxcluster_exe = ample_util.find_exe('maxcluster')


        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return
    """
    def testEnsemblingPercent(self):
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh5")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")

        num_clusters = 1
        cluster_method = 'spicker'
        percent_truncation = 5
        truncation_method = "percent"
        ensembles = ensembler.generate_ensembles(models,
                                                 cluster_method=cluster_method,
                                                 cluster_exe=self.spicker_exe,
                                                 num_clusters=num_clusters,
                                                 percent_truncation=percent_truncation,
                                                 truncation_method=truncation_method,
                                                 work_dir=work_dir)
        
        # Below tested with ccp4 6.5.010 on osx 10.9.5
        eref = ['c1_t100_r2_allatom.pdb', 'c1_t100_r2_polyAla.pdb', 'c1_t100_r2_reliable.pdb', 'c1_t100_r3_allatom.pdb',
                 'c1_t100_r3_polyAla.pdb', 'c1_t100_r3_reliable.pdb', 'c1_t19_r1_allatom.pdb', 'c1_t19_r1_polyAla.pdb',
                 'c1_t19_r1_reliable.pdb', 'c1_t24_r1_allatom.pdb', 'c1_t24_r1_polyAla.pdb', 'c1_t24_r1_reliable.pdb',
                 'c1_t29_r1_allatom.pdb', 'c1_t29_r1_polyAla.pdb', 'c1_t29_r1_reliable.pdb', 'c1_t34_r1_allatom.pdb',
                 'c1_t34_r1_polyAla.pdb', 'c1_t34_r1_reliable.pdb', 'c1_t39_r1_allatom.pdb', 'c1_t39_r1_polyAla.pdb',
                 'c1_t39_r1_reliable.pdb', 'c1_t44_r1_allatom.pdb', 'c1_t44_r1_polyAla.pdb', 'c1_t44_r1_reliable.pdb',
                 'c1_t44_r2_allatom.pdb', 'c1_t44_r2_polyAla.pdb', 'c1_t44_r2_reliable.pdb', 'c1_t49_r1_allatom.pdb',
                 'c1_t49_r1_polyAla.pdb', 'c1_t49_r1_reliable.pdb', 'c1_t49_r2_allatom.pdb', 'c1_t49_r2_polyAla.pdb',
                 'c1_t49_r2_reliable.pdb', 'c1_t54_r1_allatom.pdb', 'c1_t54_r1_polyAla.pdb', 'c1_t54_r1_reliable.pdb',
                 'c1_t54_r2_allatom.pdb', 'c1_t54_r2_polyAla.pdb', 'c1_t54_r2_reliable.pdb', 'c1_t59_r1_allatom.pdb',
                 'c1_t59_r1_polyAla.pdb', 'c1_t59_r1_reliable.pdb', 'c1_t59_r2_allatom.pdb', 'c1_t59_r2_polyAla.pdb',
                 'c1_t59_r2_reliable.pdb', 'c1_t64_r1_allatom.pdb', 'c1_t64_r1_polyAla.pdb', 'c1_t64_r1_reliable.pdb',
                 'c1_t64_r2_allatom.pdb', 'c1_t64_r2_polyAla.pdb', 'c1_t64_r2_reliable.pdb', 'c1_t69_r1_allatom.pdb',
                 'c1_t69_r1_polyAla.pdb', 'c1_t69_r1_reliable.pdb', 'c1_t69_r2_allatom.pdb', 'c1_t69_r2_polyAla.pdb',
                 'c1_t69_r2_reliable.pdb', 'c1_t75_r1_allatom.pdb', 'c1_t75_r1_polyAla.pdb', 'c1_t75_r1_reliable.pdb',
                 'c1_t75_r2_allatom.pdb', 'c1_t75_r2_polyAla.pdb', 'c1_t75_r2_reliable.pdb', 'c1_t75_r3_allatom.pdb',
                 'c1_t75_r3_polyAla.pdb', 'c1_t75_r3_reliable.pdb', 'c1_t80_r1_allatom.pdb', 'c1_t80_r1_polyAla.pdb',
                 'c1_t80_r1_reliable.pdb', 'c1_t80_r2_allatom.pdb', 'c1_t80_r2_polyAla.pdb', 'c1_t80_r2_reliable.pdb',
                 'c1_t80_r3_allatom.pdb', 'c1_t80_r3_polyAla.pdb', 'c1_t80_r3_reliable.pdb', 'c1_t85_r1_allatom.pdb',
                 'c1_t85_r1_polyAla.pdb', 'c1_t85_r1_reliable.pdb', 'c1_t85_r2_allatom.pdb', 'c1_t85_r2_polyAla.pdb',
                 'c1_t85_r2_reliable.pdb', 'c1_t85_r3_allatom.pdb', 'c1_t85_r3_polyAla.pdb', 'c1_t85_r3_reliable.pdb',
                 'c1_t90_r1_allatom.pdb', 'c1_t90_r1_polyAla.pdb', 'c1_t90_r1_reliable.pdb', 'c1_t90_r2_allatom.pdb',
                 'c1_t90_r2_polyAla.pdb', 'c1_t90_r2_reliable.pdb', 'c1_t90_r3_allatom.pdb', 'c1_t90_r3_polyAla.pdb',
                 'c1_t90_r3_reliable.pdb', 'c1_t95_r2_allatom.pdb', 'c1_t95_r2_polyAla.pdb', 'c1_t95_r2_reliable.pdb',
                 'c1_t95_r3_allatom.pdb', 'c1_t95_r3_polyAla.pdb', 'c1_t95_r3_reliable.pdb']

        self.assertEqual(sorted([os.path.basename(m) for m in ensembles]), eref)
        d = ensembler.ensembles_data[5]

        self.assertEqual(d['percent_truncation'], percent_truncation)
        self.assertEqual(d['truncation_method'], truncation_method)
        self.assertEqual(d['cluster_method'], cluster_method)
        self.assertEqual(d['num_clusters'], num_clusters)
        self.assertTrue(abs(d['truncation_variance'] - 13.035172) < 0001)
        self.assertEqual(d['ensemble_num_atoms'], 984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']), '4_S_00000002.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testEnsemblingThresh(self):
        
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler = Ensembler()

        work_dir = os.path.join(self.tests_dir, "genthresh6")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe = self.theseus_exe
        ensembler.cluster_exe = self.spicker_exe
        ensembler.subcluster_exe = self.maxcluster_exe
        
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")
        
        num_clusters = 1
        cluster_method = 'spicker'
        percent_truncation = 5
        truncation_method = "thresh"
        ensembles = ensembler.generate_ensembles(models,
                                               cluster_method=cluster_method,
                                               cluster_exe=self.spicker_exe,
                                               num_clusters=num_clusters,
                                               percent_truncation=percent_truncation,
                                               truncation_method=truncation_method,
                                               work_dir=work_dir)
        
        self.assertEqual(len(ensembles), 162, len(ensembles))
        d = ensembler.ensembles_data[5]
        
        self.assertTrue(abs(d['truncation_variance'] - 27.389253) < 0001)
        self.assertEqual(d['percent_truncation'], percent_truncation)
        self.assertEqual(d['truncation_method'], truncation_method)
        self.assertEqual(d['cluster_method'], cluster_method)
        self.assertEqual(d['num_clusters'], num_clusters)
        self.assertEqual(d['subcluster_radius_threshold'], 3)
        self.assertEqual(d['side_chain_treatment'], ALLATOM)
        self.assertEqual(d['ensemble_num_atoms'], 984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']), '5_S_00000005.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    """
    

