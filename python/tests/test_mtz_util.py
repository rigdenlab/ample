"""Test functions for python.mtz_util"""

import os
import unittest
from ample.python import mtz_util


class Test(unittest.TestCase):
    """
    Unit test
    """

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
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def test_processReflectionFile(self):
        # Get MTZ flags
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.ample_dir, "examples", "toxd-example" , "1dtx.mtz" )

        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
               'work_dir': os.getcwd(),
             }

        mtz_util.processReflectionFile( d )
        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

    def test_processReflectionFileNORFREE(self):
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
        print "GOT ",d['FREE']
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")
        
        os.unlink('uniqueify.log')
        os.unlink('2uui_sigmaa_uniqueify.mtz')
        os.unlink('2uui_sigmaa_uniqueify.log')

    def test_processReflectionFileCIF(self):
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

        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        os.unlink('cif2mtz.log')
        os.unlink('1x79-sf.mtz')
        os.unlink('mtzutils.log')
        #os.unlink('1x79-sf_dFREE.mtz')
        os.unlink('uniqueify.log')
        os.unlink('1x79-sf_uniqueify.mtz')
        os.unlink('1x79-sf_uniqueify.log')
    
    def test_processMtzLabels(self):
        # Get MTZ flags
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.testfiles_dir, "2uui_sigmaa.mtz" )
        
        FP,SIGFP,FREE=mtz_util.get_labels(mtz)
        self.assertEqual(FP,'F')
        self.assertEqual(SIGFP,'SIGF')
        self.assertEqual(FREE,None)
    
    def test_processMtzLabels2(self):
        # Get MTZ flags
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mtz = os.path.join( self.testfiles_dir, "1dtx.mtz" )
        
        FP,SIGFP,FREE=mtz_util.get_labels(mtz)
        self.assertEqual(FP,'FP')
        self.assertEqual(SIGFP,'SIGFP')
        self.assertEqual(FREE,'FreeR_flag')  

if __name__ == "__main__":
    unittest.main()