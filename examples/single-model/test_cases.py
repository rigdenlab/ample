#!/usr/bin/env ccp4-python
'''
Created on 16 Jan 2016

@author: hlfsimko
'''

import glob
import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "single-model", "input")
TEST_DICT = {}

# vanilla test
args_vanilla = [
    [ '-fasta', os.path.join(INPUT_DIR, '1ujb.fasta') ],
    [ '-mtz', os.path.join(INPUT_DIR, '1ujb-sf.mtz') ],
]

###############################################################################
#
# test ensemble creation from a single structure based on residue scores
#
###############################################################################
args_from_single_model = args_vanilla + [
    [ '-percent', '20' ],
    [ '-single_model', os.path.join(INPUT_DIR, '3c7t.pdb') ],
    [ '-truncation_scorefile', os.path.join(INPUT_DIR, '3c7t_scores.csv') ],
    [ '-truncation_scorefile_header', 'residue', 'concoord' ],                                         
]


# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_from_single_model(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 15, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['ensembles']), 0, "No ensembles produced")
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
TEST_DICT['from_single_model'] = { 'args' : args_from_single_model,
                                   'test' :  AMPLETest,
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
