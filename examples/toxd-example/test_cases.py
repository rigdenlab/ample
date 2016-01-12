#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import glob
import os
import sys

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

test_dict = {}

# vanilla test
args_vanilla = [
    [ '-fasta', 'toxd_.fasta' ],
    [ '-mtz', '1dtx.mtz' ],
    [ '-percent', '50' ]
]

###############################################################################
#
# Making models
#
###############################################################################

args_rosetta_modelling = args_vanilla + [
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-frags_3mers', 'aat000_03_05.200_v1_3' ],
    [ '-frags_9mers', 'aat000_09_05.200_v1_3' ],
    [ '-nmodels', '30']
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_rosetta_modelling(self):
        self.assertIn('models_dir', self.AMPLE_DICT)
        nmodels = len(glob.glob(self.AMPLE_DICT['models_dir']+"/*.pdb"))
        self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['ensembles']), 0, "No ensembles produced")
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
test_dict['rosetta_modelling'] = { 'args' : args_rosetta_modelling,
                                   'test' :  AMPLETest,
                                   'directory' : os.path.abspath(os.path.dirname(__file__))
                                    }

###############################################################################
#
# test from pre-existing models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_existing_models = args_vanilla + [
    [ '-models', '../../tests/testfiles/models' ],
    [ '-native_pdb', '1DTX.pdb' ],                                         
]


# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_from_existing_models(self):
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 12, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'],"Job did not succeed")
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        self.assertTrue(os.path.isfile(os.path.join(self.AMPLE_DICT['benchmark'],'results.csv')),"Missing benchmark results.csv")
        return
        
test_dict['from_existing_models'] = { 'args' : args_from_existing_models,
                                     'test' :  AMPLETest,
                                     'directory' : os.path.abspath(os.path.dirname(__file__))
                                      }

###############################################################################
#
# test from quark models 
#
###############################################################################
args_from_quark_models = args_vanilla + [
    [ '-models', '../../tests/testfiles/decoys_200.tar.gz'],
    [ '-do_mr', 'False' ]                                         
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_from_quark_models(self):
        nmodels = len(glob.glob(self.AMPLE_DICT['models_dir']+"/*.pdb"))
        self.assertEqual(nmodels, 200, "Only {0} models produced".format(nmodels))
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 18, "Incorrect number of ensembles produced: {0}".format(nensembles))
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
