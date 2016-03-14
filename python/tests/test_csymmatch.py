"""Test functions for python.csymmatch"""

import os
import unittest
from ample.python import csymmatch

class TestContacts( unittest.TestCase ):

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
        return

    def test_parse1(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        logfile = os.path.join(self.testfiles_dir, "csymmatch1.log")
        c = csymmatch.Csymmatch()
        c.parseLog(logfile=logfile, cleanup=False)
        self.assertEqual(c.changeOfHand, False )
        self.assertEqual(c.changeOfOrigin, None )
        self.assertEqual(c.chainShifts, {'a': [{'resStart': 14, 'score': 0.403697, 'resEnd': 14}, 
                                               {'resStart': 17, 'score': 0.247688, 'resEnd': 17}, 
                                               {'resStart': 32, 'score': 0.528113, 'resEnd': 44}, 
                                               {'resStart': 49, 'score': 0.268943, 'resEnd': 51}]})
        self.assertEqual(c.averageScore(), 0.36211025)
    
    def test_parse2(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        logfile = os.path.join(self.testfiles_dir, "csymmatch2.log")
        c = csymmatch.Csymmatch()
        c.parseLog(logfile=logfile, cleanup=False)
        self.assertEqual(c.changeOfHand, True)
        self.assertEqual(c.changeOfOrigin, [0, 0.5625, 0])
        self.assertEqual(c.chainShifts, {'A': [{'resStart': 1, 'score': 0.440815, 'resEnd': 27}],
                                         'B': [{'resStart': 1, 'score': 0.558538, 'resEnd': 6}]})
        self.assertEqual(c.averageScore(), 0.49967649999999997)

if __name__ == "__main__":
    unittest.main()
