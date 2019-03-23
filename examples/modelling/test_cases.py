#!/usr/bin/env ccp4-python
'''
Created on 23 March 2019

@author: jmht
'''

import os
import shutil

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

TEST_DICT = {}

###############################################################################
#
# Modelling module
#
###############################################################################

# Specify the arguments to the ensembler module to run this test case
args_default =  [
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-rosetta_flagsfile', 'flags' ],
    [ '-nmodels', '3' ]
]

class AMPLETest(AMPLEBaseTest):
    def test_modelling_default(self):
        self.assertEqual(3, len(self.AMPLE_DICT['models']))
        return
    
def copy_files(run_dir):
    """Copy inputs into run directory"""
    for f in ['toxd_.fasta', 'aat000_03_05.200_v1_3', 'aat000_09_05.200_v1_3']:
        fpath = os.path.join(SHARE_DIR, "examples", "modelling", "input", f)
        shutil.copy(fpath, os.path.join(run_duir, f))

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['modelling_default'] = {'args' : args_default,
                                  'test' :  AMPLETest,
                                  'setup' :  copy_files
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
