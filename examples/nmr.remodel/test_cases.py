#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import cPickle
import os
import sys
import unittest

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

test_dict = {}

###############################################################################
#
# NMR Remodelling
#
###############################################################################
args =  [
    [ '-name', '1t00' ],
    [ '-fasta', '1T00.fasta' ],
    [ '-mtz', '1t00.mtz' ],
    [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    [ '-nmr_model_in', '2DIZ.pdb' ],
    [ '-nmr_remodel', 'True' ],
    [ '-frags_3mers', '1t00.200.3mers' ],
    [ '-frags_9mers', '1t00.200.9mers' ],
    [ '-nmr_process', '1' ]
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_nmr_remodel(self):
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
test_dict['nmr_remodel'] = { 'args' : args,
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
