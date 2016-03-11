#!/usr/bin/env ccp4-python
'''
Created on 25 Feb 2016

@author: hlfsimko
'''

import os
import sys

#AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
#sys.path.append(os.path.join(AMPLE_DIR,'python'))
from ample.tests import test_funcs

test_dict = {}

# Universal args
args_universal = [
        [ '-fasta', '1aba.fasta' ],
        [ '-mtz', '1aba-sf.mtz' ],
]

###############################################################################
#
# Import models
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_import_models = args_universal + [
        [ '-do_mr', 'False' ],
        [ '-models', 'models' ],
        [ '-num_clusters', '3' ],
]

class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_import_models(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
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
test_dict['import_models'] = { 'args' : args_import_models,
                               'test' :  AMPLETest,
                               'directory' : os.path.abspath(os.path.dirname(__file__))
}

###############################################################################
#
# Import cluster
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_import_cluster = args_universal + [
        [ '-cluster_dir', 'models' ],
        [ '-do_mr', 'False' ],
]

class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_import_cluster(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertEqual(150, len(self.AMPLE_DICT['ensembles_data']))
        return

# Add everything to the test_dict - the key is used to name the script and run directory
test_dict['import_cluster'] = { 'args' : args_import_cluster,
                                'test' :  AMPLETest,
                                'directory' : os.path.abspath(os.path.dirname(__file__))
}

###############################################################################
#
# Import ensembles
#
###############################################################################

# Specify the arguments to AMPLE to run this test case
args_import_ensembles = args_universal + [
        [ '-ensembles', 'ensembles/' ],
        [ '-native_pdb', '1aba.pdb' ],
        [ '-shelx_cycles', '1' ],
        [ '-use_arpwarp', 'False' ],
        [ '-use_buccaneer', 'False' ],
]

class AMPLETest(test_funcs.AMPLEBaseTest):
    def test_import_ensembles(self):
        self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
        self.assertIn('ensembles', self.AMPLE_DICT)
        self.assertGreater(len(self.AMPLE_DICT['mrbump_results']), 0, "No MRBUMP results")
        self.assertTrue(self.AMPLE_DICT['success'])
        self.assertGreater(self.AMPLE_DICT['mrbump_results'][0]['SHELXE_CC'], 25,"SHELXE_CC criteria not met")
        return

# Add everything to the test_dict - the key is used to name the script and run directory
test_dict['import_ensembles'] = { 'args' : args_import_ensembles,
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
