#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import glob
import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

DIR = os.path.join(SHARE_DIR, "examples", "toxd-example")
TEST_DICT = {}
TESTFILES_DIR = os.path.join(SHARE_DIR, "testfiles")

# vanilla test
args_vanilla = [
    [ '-fasta', os.path.join(DIR, 'toxd_.fasta') ],
    [ '-mtz', os.path.join(DIR, '1dtx.mtz') ],
    [ '-percent', '50' ]
]

if not sys.platform.startswith('win'):
    ###############################################################################
    #
    # Making models
    #
    ###############################################################################

    args_rosetta_modelling = args_vanilla + [
        [ '-rosetta_dir', '/opt/rosetta-3.5' ],
        [ '-frags_3mers', os.path.join(DIR, 'aat000_03_05.200_v1_3') ],
        [ '-frags_9mers', os.path.join(DIR, 'aat000_09_05.200_v1_3') ],
        [ '-nmodels', '40']
    ]
    
    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_rosetta_modelling(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertIn('models_dir', self.AMPLE_DICT)
            nmodels = len(glob.glob(self.AMPLE_DICT['models_dir']+"/*.pdb"))
            self.assertEqual(nmodels, 40, "Only {0} models produced".format(nmodels))
            self.assertIn('ensembles', self.AMPLE_DICT)
            self.assertGreater(len(self.AMPLE_DICT['ensembles']), 0, "No ensembles produced")
            self.assertIn('mrbump_results', self.AMPLE_DICT)
            self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
            self.assertTrue(self.AMPLE_DICT['success'])
            self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
            return
            
    TEST_DICT['rosetta_modelling'] = { 'args' : args_rosetta_modelling,
                                       'test' :  AMPLETest,
                                        }

###############################################################################
#
# test from pre-existing models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_existing_models = args_vanilla + [
    [ '-models', os.path.join(TESTFILES_DIR, 'models') ],
    [ '-native_pdb', os.path.join(DIR, '1DTX.pdb') ],                                         
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_existing_models(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 12, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'],"Job did not succeed")
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        self.assertTrue(os.path.isfile(os.path.join(self.AMPLE_DICT['benchmark_dir'],'results.csv')),"Missing benchmark results.csv")
        return
        
TEST_DICT['from_existing_models'] = { 'args' : args_from_existing_models,
                                     'test' :  AMPLETest,
                                      }

###############################################################################
#
# test from quark models 
#
###############################################################################
args_from_quark_models = args_vanilla + [
    [ '-models', os.path.join(TESTFILES_DIR, 'decoys_200.tar.gz') ],
    [ '-do_mr', 'False' ]                                         
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_quark_models(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        nmodels = len(glob.glob(self.AMPLE_DICT['models_dir']+"/*.pdb"))
        self.assertEqual(nmodels, 200, "Only {0} models produced".format(nmodels))
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        ref_nensembles = 18 if self.AMPLE_DICT["use_scwrl"] else 6 # unmod ensembles without side-chains - Scwrl4 required
        self.assertEqual(nensembles, ref_nensembles, "Incorrect number of ensembles produced: {0}".format(nensembles))
        return
         
TEST_DICT['from_quark_models'] = { 'args' : args_from_quark_models,
                                   'test' :  AMPLETest,
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
