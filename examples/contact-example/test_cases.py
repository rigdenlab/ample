#!/usr/bin/env ccp4-python
'''
Created on 12 Jan 2016

@author: hlfsimko, jmht
'''

import glob
import os
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

DIR = os.path.join(SHARE_DIR, "examples", "contact-example")
TEST_DICT = {}

if not sys.platform.startswith('win'):
    # vanilla test
    args_vanilla = [
        [ '-fasta', os.path.join(DIR, 'toxd_.fasta') ],
        [ '-mtz', os.path.join(DIR, '1dtx.mtz') ],
        [ '-percent', '50' ],
        [ '-frags_3mers', os.path.join(DIR, 'aat000_03_05.200_v1_3') ],
        [ '-frags_9mers', os.path.join(DIR, 'aat000_09_05.200_v1_3') ],
        [ '-psipred_ss2', os.path.join(DIR, 'toxd_.psipred_ss2') ],
        [ '-nmodels', '30' ],
        [ '-do_mr', 'False' ],
        [ '-rosetta_dir', '/opt/rosetta-3.5' ],
    ]
    
    ###############################################################################
    #
    # Making models WITH contacts
    #
    ###############################################################################
    
    args_rosetta_contacts = args_vanilla + [
        [ '-contact_file', os.path.join(DIR, 'toxd_.pconsc2.CASPRR') ],
        [ '-energy_function', 'FADE_default' ],
    ]
    
    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_rosetta_contacts(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertEqual(self.AMPLE_DICT['energy_function'], "FADE_default")
            self.assertIn('restraints_file', self.AMPLE_DICT)
            self.assertTrue(os.path.exists(self.AMPLE_DICT['restraints_file']))
            self.assertIn('contact_map', self.AMPLE_DICT)
            self.assertTrue(os.path.isfile(self.AMPLE_DICT['contact_map']))
            self.assertIn('models_dir', self.AMPLE_DICT)
            m_dir = os.path.join(self.AMPLE_DICT['models_dir'], "..", "modelling", "job_0")
            for line in open(os.path.join(m_dir, "model_0.log"), "r"):
                if line.startswith("core.scoring.constraints.ConstraintsIO: Read in") \
                        and line.strip().endswith("constraints"):
                    nr_restraints=int(line.split()[-2])
                    break
                else: nr_restraints=-1
            self.assertGreaterEqual(nr_restraints, 0, "Restraints not read")
            self.assertEqual(nr_restraints, 59, "Different number read")
            nmodels = len(glob.glob(os.path.join(self.AMPLE_DICT['models_dir'], "*.pdb")))
            self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))
            return
     
    TEST_DICT['rosetta_contacts'] = {'args' : args_rosetta_contacts,
                                     'test' :  AMPLETest,
                                    }
    
    ###############################################################################
    #
    # Making models WITH restraints
    #
    ###############################################################################
    
    args_rosetta_restraints = args_vanilla + [
        [ '-restraints_file', os.path.join(DIR, 'toxd_.cst')]
    ]
    
    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_rosetta_contacts(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertIn('contact_map', self.AMPLE_DICT)
            self.assertTrue(os.path.isfile(self.AMPLE_DICT['contact_map']))
            self.assertIn('models_dir', self.AMPLE_DICT)
            m_dir = os.path.join(self.AMPLE_DICT['models_dir'], "..", "modelling", "job_0")
            for line in open(os.path.join(m_dir, "model_0.log"), "r"):
                if line.startswith("core.scoring.constraints.ConstraintsIO: Read in") \
                        and line.strip().endswith("constraints"):
                    nr_restraints=int(line.split()[-2])
                    break
                else: nr_restraints=-1
            self.assertGreaterEqual(nr_restraints, 0, "Restraints not read")
            self.assertEqual(nr_restraints, 49, "Different number read")
            nmodels = len(glob.glob(os.path.join(self.AMPLE_DICT['models_dir'], "*.pdb")))
            self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))
            return
    
    TEST_DICT['rosetta_restraints'] = {'args' : args_rosetta_restraints,
                                       'test' :  AMPLETest,
                                      }

if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
