
import shutil
import unittest

from ample.ensembler.homologs import *

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
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        
        cls.theseus_exe = ample_util.find_exe('theseus')
        cls.spicker_exe = ample_util.find_exe('spicker')
        cls.maxcluster_exe = ample_util.find_exe('maxcluster')


        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return
    
    def test_gesamt(self):
        gesamt_exe = "/opt/ccp4-devtools/install/bin/gesamt"
        if not ample_util.is_exe(gesamt_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "gesamt_test")
        alignment_file = align_gesamt(models, gesamt_exe=gesamt_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return
    
    def test_mustang(self):
        mustang_exe = "/opt/MUSTANG_v3.2.2/bin/mustang-3.2.1"
        if not ample_util.is_exe(mustang_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "mustang_test")
        alignment_file = align_mustang(models, mustang_exe=mustang_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return

if __name__ == "__main__":
    unittest.main()
