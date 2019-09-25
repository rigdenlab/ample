"""Test functions for ample.util.contact_util"""

__author__ = "Felix Simkovic"
__date__ = "06 Dec 2016"

import tempfile
import unittest

from ample.util import contact_util


class TestSubselectionAlgorithm(unittest.TestCase):
    def test_cutoff_1(self):
        data = [1.0, 0.6, 0.5, 0.4, 0.4, 0.3, 0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.cutoff(data)
        self.assertEqual([0, 1, 2, 3, 4, 5], keep)
        self.assertEqual([6, 7], throw)

    def test_cutoff_2(self):
        data = [1.0, 0.3, 0.2, 0.1, 0.6, 0.5, 0.4, 0.4]
        keep, throw = contact_util.SubselectionAlgorithm.cutoff(data)
        self.assertEqual([0, 1, 4, 5, 6, 7], keep)
        self.assertEqual([2, 3], throw)

    def test_cutoff_3(self):
        data = [0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.cutoff(data)
        self.assertEqual([], keep)
        self.assertEqual([0, 1], throw)

    def test_cutoff_4(self):
        data = [0.286, 0.287, 0.288]
        keep, throw = contact_util.SubselectionAlgorithm.cutoff(data)
        self.assertEqual([1, 2], keep)
        self.assertEqual([0], throw)

    def test_linear_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([4, 5, 6, 7], throw)

    def test_linear_2(self):
        data = [0.1, 0.2, 0.3, 0.45, 0.4, 0.5, 0.6, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.linear(data)
        self.assertEqual([7, 6, 5, 3], keep)
        self.assertEqual([4, 2, 1, 0], throw)

    def test_linear_3(self):
        data = [1.0, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([4, 5, 6], throw)

    def test_linear_4(self):
        data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.linear(data)
        self.assertEqual([6, 5, 4, 3], keep)
        self.assertEqual([2, 1, 0], throw)

    def test_linear_5(self):
        data = [1.0, 0.21, 0.4, 0.6, 0.5, 0.3, 0.1, 0.2]
        keep, throw = contact_util.SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 3, 4, 2], keep)
        self.assertEqual([5, 1, 7, 6], throw)

    def test_scaled_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.scaled(data)
        self.assertEqual([0, 1, 2, 3, 4, 5], keep)
        self.assertEqual([6, 7], throw)

    def test_scaled_2(self):
        data = [1.0, 1.0, 1.0, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.scaled(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)

    def test_scaled_3(self):
        data = [100.0, 1.0, 1.0, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.scaled(data)
        self.assertEqual([0], keep)
        self.assertEqual([1, 2, 3], throw)

    def test_ignore_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = contact_util.SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7], keep)
        self.assertEqual([], throw)

    def test_ignore_2(self):
        data = [1.0, 1.0, 1.0, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)

    def test_ignore_3(self):
        data = [100.0, 1.0, 1.0, 1.0]
        keep, throw = contact_util.SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)


class TestContactUtil(unittest.TestCase):

    _OPTD = {
        'contact_format': 'casprr',
        'restraints_format': 'rosetta',
        'sequence_format': 'fasta',
        'structure_format': 'pdb',
        'energy_function': 'FADE',
        'subselection_mode': None,
        'bbcontacts_file': tempfile.NamedTemporaryFile(delete=False).name,
        'contact_file': tempfile.NamedTemporaryFile(delete=False).name,
        'fasta': tempfile.NamedTemporaryFile(delete=False).name,
        'native_pdb': tempfile.NamedTemporaryFile(delete=False).name,
        'native_pdb_std': tempfile.NamedTemporaryFile(delete=False).name,
        'distance_to_neighbour': 10,
        'native_cutoff': 8,
        'restraints_factor': 1.5,
        'name': 'ample_test',
        'work_dir': '.',
    }

    def test_check_options_1(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['contact_file'] = None
            contact_util.ContactUtil.check_options(x)

    def test_check_options_2(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['contact_file'] = 'foo.bar'
            contact_util.ContactUtil.check_options(x)

    def test_check_options_3(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['bbcontacts_file'] = 'foo.bar'
            contact_util.ContactUtil.check_options(x)

    def test_check_options_4(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['contact_format'] = None
            contact_util.ContactUtil.check_options(x)

    def test_check_options_5(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['contact_format'] = 'foo'
            contact_util.ContactUtil.check_options(x)

    def test_check_options_6(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['energy_function'] = 'DEFAULT'
            x['restraints_format'] = 'rosetta'
            contact_util.ContactUtil.check_options(x)

    def test_check_options_7(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['energy_function'] = 'FADE'
            x['restraints_format'] = 'saint2'
            contact_util.ContactUtil.check_options(x)

    def test_check_options_8(self):
        x = TestContactUtil._OPTD.copy()
        x['subselect_mode'] = 'linear'
        contact_util.ContactUtil.check_options(x)

    def test_check_options_9(self):
        with self.assertRaises(ValueError):
            x = TestContactUtil._OPTD.copy()
            x['subselect_mode'] = 'test'
            contact_util.ContactUtil.check_options(x)


if __name__ == "__main__":
    unittest.main(verbosity=2)
