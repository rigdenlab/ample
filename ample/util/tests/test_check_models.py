"""Test facility for ample/util/check_models.py"""

__author__ = "Jens Thomas"

import glob
import iotbx.pdb
import os
import tempfile
import unittest

from ample import constants
from ample.util import pdb_edit
from ample.util import check_models


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath(os.path.dirname(__file__)) 
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')

    def test_check_models_nmr(self):
        pdbin = os.path.join(self.ample_share, 'examples', 'nmr-truncate', 'input', '2LC9.pdb')
        results = check_models.CheckModelsResult()
        check_models.check_models([pdbin], results)
        self.assertTrue(results.nmr)
        self.assertIsNone(results.error)
        self.assertFalse(results.homolog)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.single_structure)

    def test_check_models_single_structure(self):
        pdbin = os.path.join(self.testfiles_dir, '1K33_S_00000001.pdb')
        results = check_models.CheckModelsResult()
        check_models.check_models([pdbin], results)
        self.assertTrue(results.single_structure)
        self.assertIsNone(results.error)
        self.assertFalse(results.homolog)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.nmr)

    def test_check_models_single_structure_fail(self):
        """Should fail as has a non-protein chain present"""
        pdbin = os.path.join(self.testfiles_dir, '1K33.pdb')
        results = check_models.CheckModelsResult()
        check_models.check_models([pdbin], results)
        self.assertTrue(results.error)

    def Xtest_check_models_single_structure(self):
        pdbin = os.path.join(self.testfiles_dir, '1K33.pdb')
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        print("GOT OUTDIR ",os.path.abspath(outdir))
        results = check_models.CheckModelsResult()
        results.models_dir = outdir
        check_models.check_models([pdbin], results)
        self.assertTrue(results.nmr)
        self.assertIsNone(results.error)
        self.assertFalse(results.homolog)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.single_structure)
        print("GOT %s" % results)
        os.rmdir(outdir)


if __name__ == "__main__":
    unittest.main(verbosity=2)
