"""Test functions for parsers.tmscore_parser"""
import os
import unittest
from ample.parsers import tmscore_parser

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def test_parse(self):
        logfile = os.path.join(self.testfiles_dir, "tmscore.log")
        TM = tmscore_parser.TMscoreLogParser()
        TM.parse(logfile)

        self.assertEqual(173, TM.nrResiduesCommon)
        self.assertEqual(6.654, TM.rmsd)
        self.assertEqual(0.5512, TM.tm)
        self.assertEqual(0.3147, TM.maxsub)
        self.assertEqual(0.4292, TM.gdtts)
        self.assertEqual(0.2283, TM.gdtha)
        
if __name__ == "__main__":
    unittest.main()
