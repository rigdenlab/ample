#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import glob
import os
import sys

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

test_dict = {}

# vanilla test
args_vanilla = [
    [ '-fasta', '1ujb.fasta' ],
    [ '-mtz', '1ujb-sf.mtz' ],
]

###############################################################################
#
# test ensemble creation from a single structure based on residue scores
#
###############################################################################
args_from_single_model = args_vanilla + [
    [ '-percent', '15' ],
    [ '-single_model', '3c7t.pdb' ],
    [ '-truncation_scorefile', '3c7t_scores.csv' ],
    [ '-truncation_scorefile_header', 'residue', 'concoord' ],                                         
]


# Test class that holds the functions to test the RESULTS_PKL file that will be passed in
class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_from_single_model(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        nensembles = len(self.AMPLE_DICT['ensembles'])
        self.assertEqual(nensembles, 21, "Incorrect number of ensembles produced: {0}".format(nensembles))
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['ensembles']), 0, "No ensembles produced")
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return
        
test_dict['from_single_model'] = { 'args' : args_from_single_model,
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
