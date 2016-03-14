"""Test functions for python.shelxe"""

import os
import unittest
from ample.util import ample_util
from ample.util import shelxe

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        
        cls.shelxe_exe=ample_util.find_exe('shelxe')

        return

    def test_shelxe1BYZ(self):
        pdb=os.path.join(self.testfiles_dir,"1BYZ.pdb")
        mtz=os.path.join(self.testfiles_dir,"1BYZ-cad.mtz")
        mrPdb=os.path.join(self.testfiles_dir,"1BYZ_phaser_loc0_ALL_poly_ala_trunc_0.486615_rad_1_UNMOD.1.pdb")
        origin=shelxe.shelxe_origin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[0.326, 0.19, 0.275])
        return
    
    def test_shelxe1D7M(self):
        pdb=os.path.join(self.testfiles_dir,"1D7M.pdb")
        mtz=os.path.join(self.testfiles_dir,"1D7M-cad.mtz")
        mrPdb=os.path.join(self.testfiles_dir,"1D7M_phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_5.241154_rad_1_UNMOD.1.pdb")
        origin=shelxe.shelxe_origin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[-0.0, -0.0, 0.5])
        return
