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

DIR = os.path.join(SHARE_DIR, "examples", "nmr-truncate")
TEST_DICT = {}

###############################################################################
#
# NMR Truncation
#
###############################################################################
args =  [
    [ '-name', '102l' ],
    [ '-fasta', os.path.join(DIR, '102L.fasta') ],
    [ '-mtz', os.path.join(DIR, '102l.mtz') ],
    [ '-nmr_model_in', os.path.join(DIR, '2LC9.pdb') ],
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_nmr_truncate(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
TEST_DICT['nmr_truncate'] = { 'args' : args,
                              'test' :  AMPLETest,
                               }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
