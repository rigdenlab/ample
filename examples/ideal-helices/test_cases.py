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
# Ideal Helices
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_ideal_helices =  [
                       '-fasta', '2OVC.fasta',
                       '-mtz', '2OVC-cad.mtz',
                       '-ideal_helices', 'True',
]

# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(unittest.TestCase):
    RESULTS_PKL = None
    def test_ideal_helices(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL),"Missing pkl file: {0}".format(self.RESULTS_PKL))
        with open(self.RESULTS_PKL) as f: ad = cPickle.load(f)
        self.assertIn('mrbump_results', ad)
        self.assertGreater(len(ad['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(ad['success'])
        self.assertGreater(ad['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
test_dict['ideal_helices'] = { 'args' : args_ideal_helices,
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
