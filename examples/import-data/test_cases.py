#!/usr/bin/env ccp4-python
'''
Created on 25 Feb 2016

@author: hlfsimko
'''

import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "import-data", "input")
TEST_DICT = {}

# Universal args
args_universal = [
        [ '-fasta', os.path.join(INPUT_DIR, '1aba.fasta') ],
        [ '-mtz', os.path.join(INPUT_DIR, '1aba-sf.mtz') ],
]

###############################################################################
#
# Import models
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_import_models = args_universal + [
        [ '-models', os.path.join(INPUT_DIR, 'models') ],
        [ '-do_mr', 'False' ],
        [ '-num_clusters', '3' ],
]

class AMPLETest(AMPLEBaseTest):
    def test_import_models(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertEqual(3, self.AMPLE_DICT['num_clusters'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        ensembles_data = self.AMPLE_DICT['ensembles_data']
        self.assertEqual(70, len(ensembles_data))
        for i in xrange(1, 4):
            cluster_ensembles = [ens for ens in ensembles_data if ens['cluster_num']==i]
            cluster_num_models = cluster_ensembles[0]['cluster_num_models']
            switch = {1: (6, 25), 2: (4, 22), 3: (3, 23)}
            num_models, num_ensembles = switch[i]
            self.assertEqual(num_ensembles, len(cluster_ensembles))
            self.assertEqual(num_models, cluster_num_models)
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['import_models'] = { 'args' : args_import_models,
                               'test' :  AMPLETest,
}

###############################################################################
#
# Import cluster
#
###############################################################################

# As we are just running one cluster, we use all subcluster and sidechain treatments
# This also gives us a chance to test the argument processing for those options
args_import_cluster = args_universal + [
        [ '-subcluster_radius_thresholds', 1,2,3 ],
        [ '-side_chain_treatments', 'polyAla', 'reliable', 'allatom' ],
        [ '-cluster_dir', os.path.join(INPUT_DIR, 'models') ],
        [ '-do_mr', 'False' ],
]

class AMPLETest(AMPLEBaseTest):
    def test_import_cluster(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertEqual(150, len(self.AMPLE_DICT['ensembles_data']))
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['import_cluster'] = { 'args' : args_import_cluster,
                                'test' :  AMPLETest,
}

###############################################################################
#
# Import ensembles
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_import_ensembles = args_universal + [
        [ '-ensembles', os.path.join(INPUT_DIR, 'ensembles') ],
        [ '-native_pdb', os.path.join(INPUT_DIR, '1aba.pdb') ],
        [ '-shelx_cycles', '1' ],
        [ '-use_arpwarp', 'False' ],
        [ '-use_buccaneer', 'False' ],
]

class AMPLETest(AMPLEBaseTest):
    def test_import_ensembles(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
TEST_DICT['import_ensembles'] = { 'args' : args_import_ensembles,
                                  'test' :  AMPLETest,
}

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
     test_funcs.parse_args(TEST_DICT)
