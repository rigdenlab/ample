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

DIR = os.path.join(SHARE_DIR, "examples", "missing-domain.1k04")
TEST_DICT = {}

###############################################################################
#
# Missing Domain
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_missing_domain =  [
       [ '-fasta', os.path.join(DIR, '1k04_.fasta') ],
       [ '-mtz', os.path.join(DIR, '1k04_cad-unique.mtz') ],
       [ '-domain_all_chains_pdb', os.path.join(DIR, 'Known_40.pdb') ],
       [ '-missing_domain', 'True' ],
       [ '-frags_3mers', os.path.join(DIR, 'aa1k04_03_05.200_v1_3') ],
       [ '-frags_9mers', os.path.join(DIR, 'aa1k04_09_05.200_v1_3') ],
       [ '-rosetta_dir', '/opt/rosetta-3.5' ],
]

class AMPLETest(AMPLEBaseTest):
    def test_missing_domain(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('mrbump_results', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['missing_domain'] = { 'args' : args_missing_domain,
                               'test' :  AMPLETest,
                                }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
