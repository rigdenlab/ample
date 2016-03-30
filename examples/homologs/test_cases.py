#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import os

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "homologs", "input")
TEST_DICT = {}

###############################################################################
#
# Homologs
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_homologs =  [
                       [ '-fasta', os.path.join(INPUT_DIR, '3DCY.fasta') ],
                       [ '-mtz', os.path.join(INPUT_DIR, '3dcy-sf.mtz') ],
                       [ '-homologs', 'True' ],
                       [ '-models', INPUT_DIR ]
]

class AMPLETest(AMPLEBaseTest):
    def test_homologs(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['homologs'] = {'args' : args_homologs,
                         'test' :  AMPLETest,
                         }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
