"""Test functions for ample.util.contact_util"""

__author__ = "Felix Simkovic"
__date__ = "06 Dec 2016"

from ample.util import contact_util

import tempfile
import unittest


class Test(unittest.TestCase):

    _OPTD = {
        'contact_format': 'casprr',
        'restraints_format': 'rosetta',
        'sequence_format': 'fasta',
        'structure_format': 'pdb',
        'energy_function': 'FADE',
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

    def test_bbcontacts_format(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.bbcontacts_format = 'bbcontacts'
        self.assertEqual('bbcontacts', cutil.bbcontacts_format)
        with self.assertRaises(ValueError):
            cutil.bbcontacts_format = 'pconsc2'

    def test_contact_format(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.contact_format = 'pconsc2'
        self.assertEqual('pconsc2', cutil.contact_format)
        cutil.contact_format = 'casprr'
        self.assertEqual('casprr', cutil.contact_format)
        cutil.contact_format = 'pdb'
        self.assertEqual('pdb', cutil.contact_format)
        with self.assertRaises(ValueError):
            cutil.contact_format = 'fasta'

    def test_energy_function(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.restraint_format = 'rosetta'
        cutil.energy_function = 'FADE'
        self.assertEqual('FADE', cutil.energy_function)
        cutil.energy_function = 'SIGMOID_gremlin'
        self.assertEqual('SIGMOID_gremlin', cutil.energy_function)
        with self.assertRaises(ValueError):
            cutil.energy_function = 'DEFAULT'
        cutil.restraint_format = 'saint2'
        cutil.energy_function = 'DEFAULT'
        self.assertEqual('DEFAULT', cutil.energy_function)
        with self.assertRaises(ValueError):
            cutil.energy_function = 'FADE'

    def test_restraint_format(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.restraint_format = 'rosetta'
        self.assertEqual('rosetta', cutil.restraint_format)
        cutil.restraint_format = 'saint2'
        self.assertEqual('saint2', cutil.restraint_format)
        with self.assertRaises(ValueError):
            cutil.restraint_format = 'saint'

    def test_sequence_format(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.sequence_format = 'fasta'
        self.assertEqual('fasta', cutil.sequence_format)
        with self.assertRaises(ValueError):
            cutil.sequence_format = 'a3m'

    def test_structure_format(self):
        cutil = contact_util.ContactUtil(Test._OPTD)
        cutil.structure_format = 'pdb'
        self.assertEqual('pdb', cutil.structure_format)
        with self.assertRaises(ValueError):
            cutil.structure_format = 'mol'

    def test_check_options(self):
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['contact_file'] = None
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['contact_file'] = 'foo.bar'
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['bbcontacts_file'] = 'foo.bar'
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['contact_format'] = None
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['contact_format'] = 'foo'
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['energy_function'] = 'DEFAULT'
            x['restraints_format'] = 'rosetta'
            contact_util.ContactUtil.check_options(x)
        with self.assertRaises(ValueError):
            x = Test._OPTD.copy()
            x['energy_function'] = 'FADE'
            x['restraints_format'] = 'saint2'
            contact_util.ContactUtil.check_options(x)


if __name__ == "__main__":
    unittest.main(verbosity=2)