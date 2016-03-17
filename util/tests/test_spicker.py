
import glob
import os
import shutil
import unittest
from ample.util import ample_util
from ample.util import spicker
from ample.testing import test_funcs

@unittest.skipUnless(test_funcs.found_exe("spicker" + ample_util.EXE_EXT), "spicker exec missing")
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
        cls.ample_dir = os.sep.join(paths[ :-2 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        cls.spicker_exe = ample_util.find_exe('spicker' + ample_util.EXE_EXT)
    
    def test_spicker(self):
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")
        work_dir = os.path.join(self.tests_dir, "spicker")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        spickerer = spicker.Spickerer(spicker_exe=self.spicker_exe)
        spickerer.cluster(models, run_dir=work_dir)
        # This with spicker from ccp4 6.5.010 on osx 10.9.5
        names = sorted([os.path.basename(m) for m in spickerer.results[0].pdbs])
        ref = ['5_S_00000005.pdb', '4_S_00000005.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb',
                '4_S_00000003.pdb', '3_S_00000006.pdb', '3_S_00000004.pdb', '2_S_00000005.pdb',
                '2_S_00000001.pdb', '3_S_00000003.pdb', '1_S_00000005.pdb', '1_S_00000002.pdb',
                '1_S_00000004.pdb']
        self.assertEqual(names, sorted(ref)) # seem to get different results on osx
        self.assertEqual(len(names), len(ref)) 
        # Centroid of third cluster
        self.assertEqual(os.path.basename(spickerer.results[2].cluster_centroid), '5_S_00000006.pdb', 
                         "WARNING: Spicker might run differently on different operating systems")
        shutil.rmtree(work_dir)

if __name__ == "__main__":
    unittest.main()
