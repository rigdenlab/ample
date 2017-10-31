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

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "transmembrane.3LBW", "input")
TEST_DICT = {}

###############################################################################
#
# NMR Truncation
#
###############################################################################
args =  [
    [ '-name', '3LBW' ],
    [ '-fasta', os.path.join(INPUT_DIR, '3LBW.fasta') ],
    [ '-mtz', os.path.join(INPUT_DIR, '3LBW.sf.mtz') ],
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-frags_3mers', os.path.join(INPUT_DIR, '3LBW.200.3mers') ],
    [ '-frags_9mers', os.path.join(INPUT_DIR, '3LBW.200.9mers') ],
    [ '-transmembrane', 'True'],
    [ '-contact_file', os.path.join(INPUT_DIR, 'c5313da6-ce36-47ea-81ee-79925b2fa836.metapsicov.stage1.txt')],
    [ '-contact_format', 'metapsicov'],
    [ '-nmodels', '10' ],
    [ '-do_mr', 'False' ],
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(AMPLEBaseTest):
    def test_nmr_truncate(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles_data', self.AMPLE_DICT)
        return
        
TEST_DICT['transmembrane'] = { 'args' : args,
                              'test' :  AMPLETest,
                               }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
