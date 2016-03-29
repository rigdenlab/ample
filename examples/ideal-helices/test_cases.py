#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "ideal-helices", "input")
TEST_DICT = {}

###############################################################################
#
# Ideal Helices
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_ideal_helices =  [
        [ '-fasta', os.path.join(INPUT_DIR, '2OVC.fasta') ],
        [ '-mtz', os.path.join(INPUT_DIR, '2OVC-cad.mtz') ],
        [ '-ideal_helices', 'True' ],
]

class AMPLETest(AMPLEBaseTest):
    def test_ideal_helices(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['ideal_helices'] = {'args' : args_ideal_helices,
                              'test' :  AMPLETest,
                             }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
