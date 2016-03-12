#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import os
import sys

from ample.tests import test_funcs

test_dict = {}

###############################################################################
#
# Homologs
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_homologs =  [
                       [ '-fasta', '3DCY.fasta' ],
                       [ '-mtz', '3dcy-sf.mtz' ],
                       [ '-homologs', 'True' ],
                       [ '-models', '.' ]
]

class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_homologs(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
test_dict['homologs'] = { 'args' : args_homologs,
                         'test' :  AMPLETest,
                         'directory' : os.path.abspath(os.path.dirname(__file__))
                         }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(test_dict)
