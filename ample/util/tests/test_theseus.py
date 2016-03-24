
import glob
import os
import shutil
import tempfile
import unittest

from ample import constants
from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import theseus
from ample.testing import test_funcs

@unittest.skipUnless(test_funcs.found_exe("theseus" + ample_util.EXE_EXT), "theseus not found")
class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')
        cls.tests_dir = tempfile.gettempdir()
        cls.theseus_exe = ample_util.find_exe("theseus" + ample_util.EXE_EXT)

    def test_align_models(self):
        os.chdir(self.thisd)
        
        models = glob.glob(os.path.join(self.testfiles_dir,'models','*.pdb'))
        work_dir = os.path.join(self.tests_dir,'theseus_align')
        homologs = False
        rtheseus = theseus.Theseus(work_dir=work_dir,theseus_exe=self.theseus_exe)
        rtheseus.superpose_models(models,homologs=homologs)
        var_by_res = rtheseus.var_by_res()
        # Below with theseus 3.1.1 on osx 10.9.5
        ref = [(0, 1, 55.757593), (1, 2, 46.981238), (2, 3, 47.734236), (3, 4, 39.857326), (4, 5, 35.477433),
               (5, 6, 26.066719), (6, 7, 24.114493), (7, 8, 24.610988), (8, 9, 21.187142), (9, 10, 21.882375),
               (10, 11, 21.622263), (11, 12, 18.680601), (12, 13, 16.568074), (13, 14, 14.889583), (14, 15, 13.889769),
               (15, 16, 8.722903), (16, 17, 8.719501), (17, 18, 4.648107), (18, 19, 4.263961), (19, 20, 2.338545),
               (20, 21, 1.412784), (21, 22, 0.57754), (22, 23, 0.204917), (23, 24, 0.226518), (24, 25, 0.162323),
               (25, 26, 0.068066), (26, 27, 0.057023), (27, 28, 0.135811), (28, 29, 0.145613), (29, 30, 0.081845),
               (30, 31, 0.051059), (31, 32, 0.045182), (32, 33, 0.112322), (33, 34, 0.102072), (34, 35, 0.446003),
               (35, 36, 0.504418), (36, 37, 1.276947), (37, 38, 2.641781), (38, 39, 4.336794), (39, 40, 6.484846),
               (40, 41, 9.559536), (41, 42, 14.467942), (42, 43, 22.818975), (43, 44, 29.55385), (44, 45, 34.692256),
               (45, 46, 35.141769), (46, 47, 40.41399), (47, 48, 52.268871), (48, 49, 54.535848), (49, 50, 49.527155),
               (50, 51, 67.9861), (51, 52, 58.661069), (52, 53, 41.802971), (53, 54, 57.085415), (54, 55, 71.944127),
               (55, 56, 57.893953), (56, 57, 54.34137), (57, 58, 77.736775), (58, 59, 83.279371)]
        
        self.assertEqual([x.idx for x in var_by_res],[x[0] for x in ref])
        self.assertEqual([x.resSeq for x in var_by_res],[x[1] for x in ref])
        for i,(t,r) in enumerate(zip([x.variance for x in var_by_res], [x[2] for x in ref])):
            self.assertTrue(abs(t-r) < 0.0001,"Mismatch for: {0} {1} {2}".format(i,t,r))
            
        shutil.rmtree(work_dir)
    
    def test_align_models_homo(self):
        os.chdir(self.thisd)

        work_dir = os.path.join(self.tests_dir,'theseus_align_homo')
        if not os.path.isdir(work_dir): os.mkdir(work_dir)
        pdb_list = [ '1D7M.pdb', '1GU8.pdb', '2UUI.pdb', '1K33.pdb' ,'1BYZ.pdb' ]
        models = []
        tokeep_idx = [ i for i in range(12) ]
        for pdb in pdb_list:
            pdbin = os.path.join(self.testfiles_dir,pdb)
            name = os.path.splitext(pdb)[0]
            pdbout = os.path.join(self.testfiles_dir,"{0}_cut.pdb".format(name))
            pdb_edit.select_residues(pdbin, pdbout, tokeep_idx=tokeep_idx)
            models.append(pdbout)

        homologs = True
        rtheseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        rtheseus.superpose_models(models, homologs=homologs)
        var_by_res = rtheseus.var_by_res()
        # Below with theseus 3.1.1 on osx 10.9.5
        ref  = [(0, 243, 8.049061), (1, 244, 2.614031), (2, 245, 1.343609), (3, 246, 2.261761), (4, 247, 1.112115),
                (5, 248, 0.574936), (6, 249, 0.03114), (7, 250, 0.002894), (8, 251, 0.002314), (9, 252, 0.002174),
                (10, 253, 0.016252), (11, 254, 0.109965)]

        self.assertEqual([x.idx for x in var_by_res],[x[0] for x in ref])
        self.assertEqual([x.resSeq for x in var_by_res],[x[1] for x in ref])
        for i,(t,r) in enumerate(zip([x.variance for x in var_by_res], [x[2] for x in ref])):
            self.assertTrue(abs(t-r) < 0.0001,"Mismatch for: {0} {1} {2}".format(i,t,r))

        self.assertTrue(all([os.path.isfile(os.path.join(work_dir,m)) for m in rtheseus.aligned_models]))
        # clean up
        for m in models: os.unlink(m)
        shutil.rmtree(work_dir)

if __name__ == "__main__":
    unittest.main()
