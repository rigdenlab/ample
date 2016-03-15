
import os
import shutil
import unittest

from ample.ensembler.homologs import align_gesamt, align_mustang
from ample.testing.test_funcs import found_exe

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')

    def test_gesamt(self):
        gesamt_exe = os.path.join(os.environ['CCP4'], "bin", "gesamt")
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "gesamt_test")
        alignment_file = align_gesamt(models, gesamt_exe=gesamt_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
    
    @unittest.skipUnless(found_exe("mustang"), "mustang exec missing")
    def test_mustang(self):
        mustang_exe = ample_util.find_exe("mustang")
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "mustang_test")
        alignment_file = align_mustang(models, mustang_exe=mustang_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)

if __name__ == "__main__":
    unittest.main()
