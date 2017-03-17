#!/usr/bin/env ccp4-python

__author__ = "Felix Simkovic & Jens Thomas"
__date__ = "12 Jan 2016"

import glob
import os
import re
import sys

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.integration_util import AMPLEBaseTest

INPUT_DIR = os.path.join(SHARE_DIR, "examples", "contact-example", "input")
TEST_DICT = {}


def parse_restraints(logfile):
    re35 = re.compile('^core.io.constraints: Read in (\d+) constraints') 
    re36 = re.compile('^core.scoring.constraints.ConstraintsIO: Read in (\d+) constraints')
    re_objects = (re35, re36)
    with open(logfile) as f:
        for line in f:
            for reo in re_objects:
                m = reo.match(line)
                if m:
                    nr_restraints = int(m.group(1))
                    return nr_restraints
    return -1


if not sys.platform.startswith('win'):
    # vanilla test
    args_vanilla = [
        ['-fasta', os.path.join(INPUT_DIR, 'toxd_.fasta') ],
        ['-mtz', os.path.join(INPUT_DIR, '1dtx.mtz') ],
        ['-percent', '50' ],
        ['-frags_3mers', os.path.join(INPUT_DIR, 'aat000_03_05.200_v1_3')],
        ['-frags_9mers', os.path.join(INPUT_DIR, 'aat000_09_05.200_v1_3')],
        ['-psipred_ss2', os.path.join(INPUT_DIR, 'toxd_.psipred_ss2')],
        ['-nmodels', '30'],
        ['-do_mr', 'False'],
        ['-rosetta_dir', '/opt/rosetta-3.5'],
    ]
    
    ###############################################################################
    #
    # Making models WITH contacts
    #
    ###############################################################################
    
    args_rosetta_contacts = args_vanilla + [
        ['-contact_file', os.path.join(INPUT_DIR, 'toxd_.mat')],
        ['-contact_format', 'ccmpred'],
        ['-energy_function', 'FADE_default'],
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
            m_file = os.path.join(self.AMPLE_DICT['models_dir'], "..", "modelling", "job_0", "model_0.log")
            nr_restraints = parse_restraints(m_file)
            self.assertGreaterEqual(nr_restraints, 0, "Restraints not read")
            self.assertEqual(nr_restraints, 59, "Different number read")
            nmodels = len(self.AMPLE_DICT['models'])
            self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))

    TEST_DICT['rosetta_contacts'] = {
        'args': args_rosetta_contacts,
        'test':  AMPLETest,
    }

    ###############################################################################
    #
    # Sub-selecting models WITH contacts
    #
    ###############################################################################

    args_rosetta_contacts = args_vanilla + [
        ['-contact_file', os.path.join(INPUT_DIR, 'toxd_.mat')],
        ['-contact_format', 'ccmpred'],
        ['-energy_function', 'FADE_default'],
        ['-subselect_mode', 'linear'],
    ]

    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_rosetta_contacts_subselect(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertIn('models', self.AMPLE_DICT)
            nmodels = len(glob.glob(os.path.join(self.AMPLE_DICT['models_dir'], '*pdb')))
            self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))
            nmodels_sub = len(self.AMPLE_DICT['models'])
            self.assertEqual(nmodels_sub, 15, "Only {0} models produced".format(nmodels))

    TEST_DICT['rosetta_contacts_subselect'] = {
        'args': args_rosetta_contacts,
        'test':  AMPLETest,
    }
    
    ###############################################################################
    #
    # Making models WITH restraints
    #
    ###############################################################################
    
    args_rosetta_restraints = args_vanilla + [
        ['-restraints_file', os.path.join(INPUT_DIR, 'toxd_.cst')]
    ]
    
    # Test class that holds the functions to test the RESULTS_PKL file that will be passed in
    class AMPLETest(AMPLEBaseTest):
        def test_rosetta_restraints(self):
            self.assertTrue(self.AMPLE_DICT['AMPLE_finished'])
            self.assertIn('models_dir', self.AMPLE_DICT)
            m_file = os.path.join(self.AMPLE_DICT['models_dir'], "..", "modelling", "job_0", "model_0.log")
            nr_restraints = parse_restraints(m_file)
            self.assertGreaterEqual(nr_restraints, 0, "Restraints not read")
            self.assertEqual(nr_restraints, 49, "Different number read")
            nmodels = len(self.AMPLE_DICT['models'])
            self.assertEqual(nmodels, 30, "Only {0} models produced".format(nmodels))

    TEST_DICT['rosetta_restraints'] = {
        'args': args_rosetta_restraints,
        'test':  AMPLETest,
    }


if __name__ == '__main__':
    test_funcs.parse_args(TEST_DICT)
