#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import os

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

TEST_DICT = {}

###############################################################################
#
# Ensembler module
#
###############################################################################

# Specify the arguments to the ensembler module to run this test case
args_default =  [
                       [ '-models', os.path.join(SHARE_DIR, "examples", "import-data", "input", "models") ],
                       [ '-num_clusters', '3' ]
]

class AMPLETest(AMPLEBaseTest):
    def test_ensembler_default(self):
        self.assertEqual(3, self.AMPLE_DICT['num_clusters'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        ensembles_data = self.AMPLE_DICT['ensembles_data']
        self.assertEqual(243, len(ensembles_data))
        for i in xrange(1, 4):
            cluster_ensembles = [ens for ens in ensembles_data if ens['cluster_num']==i]
            cluster_num_models = cluster_ensembles[0]['cluster_num_models']
            switch = {1: (6, 93), 2: (4, 72), 3: (3, 78)}
            num_models, num_ensembles = switch[i]
            self.assertEqual(num_ensembles, len(cluster_ensembles))
            self.assertEqual(num_models, cluster_num_models)
        return

# Add everything to the test_dict - the key is used to name the script and run directory
# The ensembler integration tests are special as each case key must begin with the word 'ensembler'
TEST_DICT['ensembler_default'] = {'args' : args_default,
                                  'test' :  AMPLETest,
                                 }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
