#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import cPickle
import glob
import os
import sys
import unittest

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

test_dict = {}

# vanilla test
args_vanilla = [
    '-fasta', 'toxd_.fasta',
    '-mtz', '1dtx.mtz',
    '-percent', '50',
]

###############################################################################
#
# Making models
#
###############################################################################

args_rosetta_modelling = args_vanilla + [
    '-rosetta_dir', '/opt/rosetta-3.5',
    '-frags_3mers', 'aat000_03_05.200_v1_3',
    '-frags_9mers', 'aat000_09_05.200_v1_3',
    '-nmodels', '30',
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(unittest.TestCase):
    RESULTS_PKL = None
    def test_rosetta_modelling(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL),"Missing pkl file: {0}".format(self.RESULTS_PKL))
        with open(self.RESULTS_PKL) as f: ad = cPickle.load(f)
        self.assertIn('models_dir', ad)
        nmodels = len(glob.glob(ad['models_dir']+"/*.pdb"))
        self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))
        self.assertIn('ensembles', ad)
        self.assertGreater(len(ad['ensembles']), 0, "No ensembles produced")
        self.assertIn('mrbump_results', ad)
        self.assertGreater(len(ad['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(ad['success'])
        self.assertGreater(ad['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
test_dict['rosetta_modelling'] = { 'args' : args_rosetta_modelling,
                                   'test' :  AMPLETest,
                                   'directory' : os.path.abspath(os.path.dirname(__file__))
                                    }

###############################################################################
#
# test from pre-existing models
#
###############################################################################
args_from_existing_models = args_vanilla + [
    '-models', '../../tests/testfiles/models',
]


# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(unittest.TestCase):
    RESULTS_PKL = None
    def test_from_existing_models(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL),"Missing pkl file: {0}".format(self.RESULTS_PKL))
        with open(self.RESULTS_PKL) as f: ad = cPickle.load(f)
        self.assertIn('ensembles', ad)
        nensembles = len(ad['ensembles'])
        self.assertEqual(nensembles, 12, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('mrbump_results', ad)
        self.assertGreater(len(ad['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(ad['success'])
        self.assertGreater(ad['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
test_dict['from_existing_models'] = { 'args' : args_from_existing_models,
                                     'test' :  AMPLETest,
                                     'directory' : os.path.abspath(os.path.dirname(__file__))
                                      }

###############################################################################
#
# test from quark models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_quark_models = args_vanilla + [
    '-models', '../../tests/testfiles/decoys_200.tar.gz',
    '-native_pdb', '1DTX.pdb'                                         
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(unittest.TestCase):
    RESULTS_PKL = None
    def test_from_quark_models(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL),"Missing pkl file: {0}".format(self.RESULTS_PKL))
        with open(self.RESULTS_PKL) as f: ad = cPickle.load(f)
        self.assertIn('ensembles', ad)
        nensembles = len(ad['ensembles'])
        self.assertEqual(nensembles, 18, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('mrbump_results', ad)
        self.assertGreater(len(ad['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(ad['success'])
        self.assertGreater(ad['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
         
test_dict['from_quark_models'] = { 'args' : args_from_quark_models,
                                   'test' :  AMPLETest,
                                   'directory' : os.path.abspath(os.path.dirname(__file__))
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(test_dict)
