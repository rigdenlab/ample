"""Test functions for python.shelxe"""

import os
import unittest
from ample import constants
from ample.util import ample_util
from ample.util import shelxe
from ample.testing import test_funcs

@unittest.skipUnless(test_funcs.found_exe("shelxe" + ample_util.EXE_EXT),  "shelxe exec missing")
class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')

    def test_shelxe_1BYZ(self):
        shelxe_exe = ample_util.find_exe('shelxe' + ample_util.EXE_EXT)
        pdb = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        mtz = os.path.join(self.testfiles_dir, "1BYZ-cad.mtz")
        mrPdb = os.path.join(self.testfiles_dir, 
                            "1BYZ_phaser_loc0_ALL_poly_ala_trunc_0.486615_rad_1_UNMOD.1.pdb")
        mrinfo = shelxe.MRinfo(shelxe_exe, pdb, mtz)
        mrinfo.analyse(mrPdb)
        
        self.assertEqual(mrinfo.originShift,[0.326, 0.19, 0.275])
        # Need MRBUMP to be updated
        #self.assertEqual(mrinfo.MPE,78.5)
        self.assertEqual(mrinfo.wMPE,74.5)
        
        os.unlink('shelxe-input.hkl')
        os.unlink('shelxe-input.ent')
    
    def test_shelxe_1D7M(self):
        shelxe_exe = ample_util.find_exe('shelxe' + ample_util.EXE_EXT)
        pdb = os.path.join(self.testfiles_dir, "1D7M.pdb")
        mtz = os.path.join(self.testfiles_dir, "1D7M-cad.mtz")
        mrPdb = os.path.join(self.testfiles_dir, 
                            "1D7M_phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_5.241154_rad_1_UNMOD.1.pdb")
        
        mrinfo = shelxe.MRinfo(shelxe_exe, pdb, mtz)
        mrinfo.analyse(mrPdb)
        
        self.assertEqual(mrinfo.originShift,[-0.0, -0.0, 0.5])
        # Need MRBUMP to be updated
        #self.assertEqual(mrinfo.MPE,65.5)
        self.assertEqual(mrinfo.wMPE,57.4)
        
        os.unlink('shelxe-input.hkl')
        os.unlink('shelxe-input.ent')
             
if __name__ == "__main__":
    unittest.main()
