#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import os
import sys
import unittest

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "nmr-remodel", "input")
TEST_DICT = {}

if not sys.platform.startswith('win'):
    ###############################################################################
    #
    # NMR Remodelling
    #
    ###############################################################################
    args =  [
        [ '-name', '1t00' ],
        [ '-fasta', os.path.join(INPUT_DIR, '1T00.fasta') ],
        [ '-mtz', os.path.join(INPUT_DIR,'1t00.mtz') ],
        [ '-rosetta_dir', '/opt/rosetta-3.5' ],
        [ '-nmr_model_in', os.path.join(INPUT_DIR,'2DIZ.pdb') ],
        [ '-nmr_remodel', 'True' ],
        [ '-percent_fixed_intervals', 20 ],
        [ '-frags_3mers', os.path.join(INPUT_DIR, '1t00.200.3mers') ],
        [ '-frags_9mers', os.path.join(INPUT_DIR, '1t00.200.9mers') ],
        [ '-nmr_process', '1' ]
    ]
    
    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_nmr_remodel(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertIn('mrbump_results', self.AMPLE_DICT)
            self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
            self.assertTrue(self.AMPLE_DICT['success'])
            self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
            return
            
    TEST_DICT['nmr_remodel'] = { 'args' : args,
                                 'test' :  AMPLETest,
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
