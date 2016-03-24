
import os
import shutil
import unittest

from ample.ensembler.homologs import align_gesamt, align_mustang
from ample.util import ample_util
from ample.testing import constants
from ample.testing import test_funcs

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_dir = constants.AMPLE_DIR
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')

    @unittest.skipUnless(test_funcs.found_exe("gesamt" + ample_util.EXE_EXT), "gesamt exec missing")
    def test_gesamt(self):
        gesamt_exe = ample_util.find_exe("gesamt" + ample_util.EXE_EXT)
        pdb_list = [ '1ujbA.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "gesamt_test")
        alignment_file = align_gesamt(models, gesamt_exe=gesamt_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
    
    @unittest.skipUnless(test_funcs.found_exe("mustang" + ample_util.EXE_EXT), "mustang exec missing")
    def test_mustang(self):
        mustang_exe = ample_util.find_exe("mustang" + ample_util.EXE_EXT)
        pdb_list = [ '1ujbA.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "mustang_test")
        alignment_file = align_mustang(models, mustang_exe=mustang_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)

if __name__ == "__main__":
    unittest.main()
