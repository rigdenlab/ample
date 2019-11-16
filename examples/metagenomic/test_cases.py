#!/usr/bin/env ccp4-python

__author__ = "Adam Simpkin"
__date__ = "25 Oct 2019"
__version__ = "0.1"

import glob
import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "metagenomic_example", "input")
TEST_DICT = {}
TESTFILES_DIR = os.path.join(SHARE_DIR, "testfiles")

# vanilla test
args_vanilla = [
    [ '-fasta', os.path.join(INPUT_DIR, '5edl.fasta') ],
    [ '-mtz', os.path.join(INPUT_DIR, '5edl.mtz') ],
]


###############################################################################
#
# test from pre-existing models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_existing_models = args_vanilla + [
    [ '-models', os.path.join(TESTFILES_DIR, 'metagenomic_models') ],
    [ '-native_pdb', os.path.join(INPUT_DIR, '5edl.pdb') ],                  
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_existing_models(self):

        nclusters = 10
        radii = set([1, 3])
        sidechains = set(['polyala'])
        self.assertEqual(nclusters, self.AMPLE_DICT['num_clusters'])
        self.assertEqual(list(radii), self.AMPLE_DICT['subcluster_radius_thresholds'])
        self.assertEqual(list(sidechains), self.AMPLE_DICT['side_chain_treatments'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        nensembles = len(self.AMPLE_DICT['ensembles'])
        # Was 4 changed 19/10/17 on linux CCP4 7.0.44
        # Was 3 changed 10/5/18 on linux CCP4 7.0.55
        #self.assertEqual(nensembles, 4, "Incorrect number of ensembles produced: {0}".format(nensembles))
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
        SHELXE_CC = 15 # use lower as we can't guarantee it will succeed but should check something sensible is produced
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], SHELXE_CC,"SHELXE_CC criteria not met")
        # Check some random benchmarking data
        self.assertTrue(os.path.isfile(os.path.join(self.AMPLE_DICT['benchmark_dir'],'results.csv')),"Missing benchmark results.csv")
        result = self.AMPLE_DICT['benchmark_results'][0]
        self.assertTrue(result['MR_program'], 'PHASER')
        #self.assertTrue(result['native_pdb_solvent_content'], 47.67)

TEST_DICT['from_existing_models'] = {
    'args': args_from_existing_models,
    'test':  AMPLETest,
}

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
