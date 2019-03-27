#!/usr/bin/env ccp4-python

__author__ = "Jens Thomas"
__date__ = "29 Dec 2015"
__version__ = "0.1"

import glob
import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "toxd-example", "input")
TEST_DICT = {}
TESTFILES_DIR = os.path.join(SHARE_DIR, "testfiles")

# vanilla test
args_vanilla = [
    [ '-fasta', os.path.join(INPUT_DIR, 'toxd_.fasta') ],
    [ '-mtz', os.path.join(INPUT_DIR, '1dtx.mtz') ],
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
        [ '-frags_3mers', os.path.join(INPUT_DIR, 'aat000_03_05.200_v1_3') ],
        [ '-frags_9mers', os.path.join(INPUT_DIR, 'aat000_09_05.200_v1_3') ],
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
            # Can't guarantee that this case will succeed as it fails ~ 40% of the time
            #self.assertTrue(self.AMPLE_DICT['success'])
            self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")

    TEST_DICT['rosetta_modelling'] = {
        'args': args_rosetta_modelling,
        'test':  AMPLETest,
    }

###############################################################################
#
# test from pre-existing models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_existing_models = args_vanilla + [
    [ '-models', os.path.join(TESTFILES_DIR, 'models') ],
    [ '-native_pdb', os.path.join(INPUT_DIR, '1DTX.pdb') ],                  
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_existing_models(self):

        nclusters = 10
        radii = set([1, 3])
        sidechains = set(['polyala'])
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertEqual(nclusters, self.AMPLE_DICT['num_clusters'])
        self.assertEqual(list(radii), self.AMPLE_DICT['subcluster_radius_thresholds'])
        self.assertEqual(list(sidechains), self.AMPLE_DICT['side_chain_treatments'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        # Was 4 changed 19/10/17 on linux CCP4 7.0.44
        # Was 3 changed 10/5/18 on linux CCP4 7.0.55
        self.assertEqual(nensembles, 4, "Incorrect number of ensembles produced: {0}".format(nensembles))
        ensembles_data = self.AMPLE_DICT['ensembles_data']
        # We need to make sure that only the set sidechains and radii are present in the ensembles, so we need to check that what we
        # find is a subset of the allowed values
        for i in xrange(1, len(ensembles_data) + 1):
            cluster_ensembles = [ens for ens in ensembles_data if ens['cluster_num']==i]
            side_chain_treatments = set([ens['side_chain_treatment'] for ens in ensembles_data if ens['cluster_num']==i ]  )
            subcluster_radius_thresholds = set([ens['subcluster_radius_threshold'] for ens in ensembles_data if ens['cluster_num']==i ]  )
            self.assertTrue(side_chain_treatments.issubset(sidechains), side_chain_treatments)
            self.assertTrue(subcluster_radius_thresholds.issubset(radii), subcluster_radius_thresholds)
        
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'],"Job did not succeed")
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        # Check some random benchmarking data
        self.assertTrue(os.path.isfile(os.path.join(self.AMPLE_DICT['benchmark_dir'],'results.csv')),"Missing benchmark results.csv")
        result = self.AMPLE_DICT['benchmark_results'][0]
        self.assertTrue(result['MR_program'], 'PHASER')
        self.assertTrue(result['native_pdb_solvent_content'], 47.67)
        self.assertGreater(result['PHASER_TFZ'], 5.0, "Low PHASER_TFZ: {0}".format(result['PHASER_TFZ']))

TEST_DICT['from_existing_models'] = {
    'args': args_from_existing_models,
    'test':  AMPLETest,
}

###############################################################################
#
# test from pre-existing models with classic_mode- to check we can still run in
# the old modality
#
###############################################################################
args_from_existing_models_classic = args_vanilla + [
                                                    [ '-models', os.path.join(TESTFILES_DIR, 'models') ],
                                                    [ '-classic_mode', True ],
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_existing_models_classic(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 12, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'],"Job did not succeed")
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")

TEST_DICT['from_existing_models_classic'] = {
    'args': args_from_existing_models_classic,
    'test':  AMPLETest,
}

###############################################################################
#
# test from quark models 
#
###############################################################################
args_from_quark_models = args_vanilla + [
    [ '-models', os.path.join(TESTFILES_DIR, 'decoys.tar.gz') ],
    [ '-do_mr', 'False' ]                                         
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_quark_models(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        nmodels = len(glob.glob(
            os.path.join(self.AMPLE_DICT['models_dir'], "*.pdb")
        ))
        self.assertEqual(nmodels, 200, "Only {0} models produced".format(nmodels))
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        # number was 35 - changed 19/10/17 with CCP4 7.0.44 linux
        # number was 30 - changed 10/5/18 with CCP4 7.0.55 linux
        ref_nensembles = 18 if self.AMPLE_DICT["use_scwrl"] else 35 # unmod ensembles without side-chains - Scwrl4 required
        msg = "Incorrect number of ensembles produced: {0} vs {1}".format(nensembles, ref_nensembles)
        self.assertEqual(nensembles, ref_nensembles, msg) 

TEST_DICT['from_quark_models'] = {
    'args': args_from_quark_models,
    'test':  AMPLETest,
}

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
