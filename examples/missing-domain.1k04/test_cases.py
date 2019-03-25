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

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "missing-domain.1k04", "input")
TEST_DICT = {}

if not sys.platform.startswith('win'):
    ############################################################################
    #
    # Missing Domain
    #
    ############################################################################
    
    # Specify the arguments to AMPLE to run this test case
    args_missing_domain =  [
           [ '-fasta', os.path.join(INPUT_DIR, '1k04_.fasta') ],
           [ '-mtz', os.path.join(INPUT_DIR, '1k04_cad-unique.mtz') ],
           [ '-domain_all_chains_pdb', os.path.join(INPUT_DIR, 'Known_40.pdb') ],
           [ '-missing_domain', 'True' ],
           [ '-frags_3mers', os.path.join(INPUT_DIR, 'aa1k04_03_05.200_v1_3') ],
           [ '-frags_9mers', os.path.join(INPUT_DIR, 'aa1k04_09_05.200_v1_3') ],
           [ '-rosetta_dir', '/opt/rosetta-3.5' ],
           [ '-nmodels', '100'],
           [ '-do_mr', 'False']
    ]
    
    class AMPLETest(AMPLEBaseTest):
        def test_missing_domain(self):
            self.assertEquals(len(self.AMPLE_DICT['models']), 100)
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
