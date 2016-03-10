"""Test functions for ensembler.subcluster_util.py"""

import unittest

from ample.ensembler import subcluster_util

class Test(unittest.TestCase):
        
    def testSliceSubcluster(self):
        '''Written by hlasimpk'''
        
        cluster_files = ["model_1.pdb", "model_2.pdb", "model_3.pdb",
                         "model_4.pdb", "model_5.pdb", "model_6.pdb",
                         "model_7.pdb", "model_8.pdb", "model_9.pdb",
                         "model_10.pdb", "model_11.pdb", "model_12.pdb",
                         "model_13.pdb", "model_14.pdb", "model_15.pdb",
                         "model_16.pdb", "model_17.pdb", "model_18.pdb",
                         "model_19.pdb", "model_20.pdb"]
        
        ########################################################################
        # Test Case 1
        pdbs = slice_subcluster(cluster_files, [], 30, 1, [1, 2, 3])        
        ref_pdbs = ["model_1.pdb", "model_2.pdb", "model_3.pdb", "model_4.pdb", 
                    "model_5.pdb", "model_6.pdb", "model_7.pdb", "model_8.pdb", 
                    "model_9.pdb", "model_10.pdb", "model_11.pdb", "model_12.pdb", 
                    "model_13.pdb", "model_14.pdb", "model_15.pdb", "model_16.pdb", 
                    "model_17.pdb", "model_18.pdb", "model_19.pdb", "model_20.pdb"]
        self.assertEqual(ref_pdbs, pdbs)
        
        ########################################################################
        # Test Case 2
        pdbs = slice_subcluster(cluster_files, [], 15, 1, [1, 2, 3])
        ref_pdbs = ["model_1.pdb", "model_2.pdb", "model_3.pdb", "model_4.pdb", 
                    "model_5.pdb", "model_6.pdb", "model_7.pdb", "model_8.pdb", 
                    "model_9.pdb", "model_10.pdb", "model_11.pdb", "model_12.pdb", 
                    "model_13.pdb", "model_14.pdb", "model_15.pdb"]
        self.assertEqual(ref_pdbs, pdbs)
        
        ########################################################################
        # Test Case 3
        previous_clusters = [["model_1.pdb", "model_2.pdb", "model_3.pdb",
                              "model_4.pdb", "model_5.pdb", "model_6.pdb",
                              "model_7.pdb", "model_8.pdb", "model_9.pdb",
                              "model_10.pdb", "model_11.pdb", "model_12.pdb",
                              "model_13.pdb", "model_14.pdb", "model_15.pdb"]]
        pdbs = slice_subcluster(cluster_files, previous_clusters, 
                                15, 3, [1, 2, 3])
        ref_pdbs = ["model_6.pdb", "model_7.pdb", "model_8.pdb", "model_9.pdb", 
                    "model_10.pdb", "model_11.pdb", "model_12.pdb", "model_13.pdb", 
                    "model_14.pdb", "model_15.pdb", "model_16.pdb", "model_17.pdb", 
                    "model_18.pdb", "model_19.pdb", "model_20.pdb"]
        self.assertEqual(ref_pdbs, pdbs)
        
        ########################################################################
        # Test Case 4
        previous_clusters = [["model_1.pdb", "model_2.pdb", "model_3.pdb",
                              "model_4.pdb", "model_5.pdb", "model_6.pdb",
                              "model_7.pdb", "model_8.pdb", "model_9.pdb",
                              "model_10.pdb", "model_11.pdb", "model_12.pdb",
                              "model_13.pdb", "model_14.pdb", "model_15.pdb",]]
        pdbs = slice_subcluster(cluster_files, previous_clusters,
                                15, 2, [1, 2, 3, 4, 5])
        ref_pdbs = ["model_2.pdb", "model_3.pdb", "model_4.pdb", "model_5.pdb", 
                    "model_6.pdb", "model_7.pdb", "model_8.pdb", "model_9.pdb", 
                    "model_10.pdb", "model_11.pdb", "model_12.pdb", "model_13.pdb",
                    "model_14.pdb", "model_15.pdb", "model_16.pdb"]
        self.assertEqual(ref_pdbs, pdbs)
    
    def testPickNmodels(self):
        '''Written by hlasimpk'''
        
        models = ['model_3.pdb', 'model_16.pdb', 'model_15.pdb',
                  'model_11.pdb', 'model_5.pdb', 'model_4.pdb',
                  'model_18.pdb', 'model_2.pdb', 'model_10.pdb',
                  'model_19.pdb', 'model_7.pdb', 'model_8.pdb',
                  'model_14.pdb', 'model_12.pdb', 'model_9.pdb',
                  'model_13.pdb', 'model_20.pdb', 'model_1.pdb',
                  'model_17.pdb', 'model_6.pdb']
        
        ########################################################################
        # Test Case 1   
        pdbs = pick_nmodels(models, [], 30)
        ref_pdbs = ('model_1.pdb', 'model_10.pdb', 'model_11.pdb', 'model_12.pdb', 
                    'model_13.pdb', 'model_14.pdb', 'model_15.pdb', 'model_16.pdb', 
                    'model_17.pdb', 'model_18.pdb', 'model_19.pdb', 'model_2.pdb',
                    'model_20.pdb', 'model_3.pdb', 'model_4.pdb', 'model_5.pdb', 
                    'model_6.pdb', 'model_7.pdb', 'model_8.pdb', 'model_9.pdb')
        self.assertEqual(ref_pdbs, pdbs)
        
        ########################################################################
        # Test Case 2
        clusters = [('model_1.pdb', 'model_10.pdb', 'model_11.pdb', 'model_12.pdb', 
                     'model_13.pdb', 'model_14.pdb', 'model_15.pdb', 'model_16.pdb', 
                     'model_17.pdb', 'model_18.pdb', 'model_19.pdb', 'model_2.pdb',
                     'model_20.pdb', 'model_3.pdb', 'model_4.pdb', 'model_5.pdb', 
                     'model_6.pdb', 'model_7.pdb', 'model_8.pdb', 'model_9.pdb')]
        pdbs = pick_nmodels(models, clusters, 30)
        ref_pdbs = None
        self.assertEqual(ref_pdbs, pdbs)