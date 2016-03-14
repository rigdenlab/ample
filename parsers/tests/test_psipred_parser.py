"""Test functions for parsers.psipred_parser"""
import os
import unittest
from ample.parsers import psipred_parser

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd = os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def test_parse(self):
        ss2file = os.path.join(self.testfiles_dir, "1aba_.psipred_ss2")
        PA = psipred_parser.PsipredSs2Parser(ss2file)
        ref_ss2 = "CEEEEEECCCCCCCCHHHHHHHHHHHCCCCEEEEEECCCCCCCCHHHHHHHHHHHCCCCCCCCCCCEEEEECCEEEECHHHHHHHHC"
        ss2 = "".join([i.ss for i in PA.residues])
        self.assertEqual(ref_ss2, ss2 )
        
if __name__ == "__main__":
    unittest.main()
