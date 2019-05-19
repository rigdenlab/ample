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

##################################################################################
#
# Modelling module
#
# REM: all test cases must begin with the name modelling to be run by the module
#
##################################################################################

#
# Run an arbitrary flagsfile
#
test_name = 'modelling_default'
args_default =  [
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-rosetta_flagsfile', os.path.join(SHARE_DIR, 'examples', 'modelling', 'input', 'flags')],
    [ '-nmodels', '3' ]
]

class AMPLETest1(TestCase):
    def test_modelling_default(self):
        models_dir = os.path.join(os.getcwd(), test_name, "models")
        nmodels = len(glob.glob(os.path.join(models_dir, "*.pdb")))
        self.assertEqual(nmodels, 3, "Only {} models in directory {}".format(nmodels, models_dir))
        return
    
def copy_files(run_dir):
    """Copy inputs into run directory"""
    for f in ['toxd_.fasta', 'aat000_03_05.200_v1_3', 'aat000_09_05.200_v1_3']:
        fpath = os.path.join(SHARE_DIR, "examples", "toxd-example", "input", f)
        shutil.copy(fpath, os.path.join(run_dir, f))

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT[test_name] = {'args' : args_default,
                        'test' :  AMPLETest1,
                        'setup' :  copy_files
                       }
#
# Multimer Modelling
#
test_name = 'modelling_multimer'
args_multimer =  [
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-fasta', os.path.join(SHARE_DIR, 'examples', 'toxd-example', 'input', 'toxd_.fasta')],
    [ '-frags_3mers', os.path.join(SHARE_DIR, 'examples', 'toxd-example', 'input', 'aat000_03_05.200_v1_3')],
    [ '-frags_9mers', os.path.join(SHARE_DIR, 'examples', 'toxd-example', 'input', 'aat000_09_05.200_v1_3')],
    [ '-multimer_modelling', 'dimer'],
    [ '-nmodels', '1' ]
]

class AMPLETest2(TestCase):
    def test_modelling_multimer(self):
        models_dir = os.path.join(os.getcwd(), test_name, "models")
        models = glob.glob(os.path.join(models_dir, "*.pdb"))
        self.assertEqual(len(models), 1)
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT[test_name] = {'args' : args_multimer,
                        'test' :  AMPLETest2,
                        'setup' :  copy_files
                       }

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
