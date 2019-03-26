#!/usr/bin/env ccp4-python
'''
Created on 23 March 2019

@author: jmht
'''

import glob
import os
import shutil
from unittest import TestCase

from ample.constants import SHARE_DIR
from ample.testing import test_funcs

TEST_DICT = {}

###############################################################################
#
# Modelling module
#
###############################################################################

# Specify the arguments to the ensembler module to run this test case
test_name = 'modelling_default'
args_default =  [
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-rosetta_flagsfile', os.path.join(SHARE_DIR, 'examples', 'modelling', 'input', 'flags')],
    [ '-nmodels', '3' ]
]

class AMPLETest(TestCase):
    def test_modelling_default(self):
        models_dir = os.path.join(os.getcwd(), test_name, "models")
        models = glob.glob(os.path.join(models_dir, "*.pdb"))
        self.assertEqual(len(models), 3)
        return
    
def copy_files(run_dir):
    """Copy inputs into run directory"""
    for f in ['toxd_.fasta', 'aat000_03_05.200_v1_3', 'aat000_09_05.200_v1_3']:
        fpath = os.path.join(SHARE_DIR, "examples", "toxd-example", "input", f)
        shutil.copy(fpath, os.path.join(run_dir, f))

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT[test_name] = {'args' : args_default,
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
