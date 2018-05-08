"""Test functions for util.mtz_util"""

import os
import sys
import unittest
from ample import constants
from ample.util import mtz_util

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')
        return

    def test_process_reflection_file(self):
        # Get MTZ flags
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.ample_share, "examples", "toxd-example" , "input", "1dtx.mtz" )

        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
               'work_dir': os.getcwd(),
             }

        mtz_util.processReflectionFile( d )
        self.assertEqual('FP', d['F'], "Correct F")
        self.assertEqual('SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual('FreeR_flag', d['FREE'], "Correct FREE")

    @unittest.skipIf(sys.platform.startswith('win'), "uniqueify not executable")
    def test_process_reflection_file_NORFREE(self):
        # Get MTZ flags

        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.testfiles_dir, "2uui_sigmaa.mtz" )

        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
              'work_dir': os.getcwd(),
             }

        mtz_util.processReflectionFile( d )

        self.assertEqual( 'F', d['F'], "Correct F")
        self.assertEqual( 'SIGF', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        os.unlink('uniqueify.log')
        os.unlink('2uui_sigmaa_uniqueify.mtz')
        os.unlink('2uui_sigmaa_uniqueify.log')

    @unittest.skipIf(sys.platform.startswith('win'), "uniqueify not executable")
    def test_process_reflection_file_CIF(self):
        # Get MTZ flags
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        cif = os.path.join( self.testfiles_dir, "1x79-sf.cif" )

        d = { 'mtz'     : None,
              'sf_cif'  : cif,
              'F'       : None,
              'SIGF'    : None,
              'FREE'    : None,
              'work_dir': os.getcwd(),
             }

        mtz_util.processReflectionFile( d )

        self.assertEqual('FP', d['F'], "Correct F")
        self.assertEqual('SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual('FreeR_flag', d['FREE'], "Correct FREE")

        os.unlink('cif2mtz.log')
        os.unlink('1x79-sf.mtz')
        os.unlink('1x79-sf_uniqueify.mtz')
        os.unlink('1x79-sf_uniqueify.log')

    def test_process_mtz_labels(self):
        # Get MTZ flags

        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.testfiles_dir, "2uui_sigmaa.mtz" )

        FP,SIGFP,FREE=mtz_util.get_labels(mtz)
        self.assertEqual(FP,'F')
        self.assertEqual(SIGFP,'SIGF')
        self.assertEqual(FREE,None)

    def test_process_mtz_labels2(self):
        # Get MTZ flags

        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.testfiles_dir, "1dtx.mtz" )

        FP,SIGFP,FREE=mtz_util.get_labels(mtz)
        self.assertEqual(FP,'FP')
        self.assertEqual(SIGFP,'SIGFP')
        self.assertEqual(FREE,'FreeR_flag')

    def test_max_min_resolution(self):
        mtz = os.path.join( self.ample_share, "examples", "toxd-example" , "input", "1dtx.mtz" )
        maxr, minr = mtz_util.max_min_resolution( mtz )
        self.assertAlmostEqual(maxr, 36.765, 3, "Wrong max resolution: {0}".format(maxr))
        self.assertAlmostEqual(minr, 2.201, 3, "Wrong min resolution: {0}".format(minr))

if __name__ == "__main__":
    unittest.main()
