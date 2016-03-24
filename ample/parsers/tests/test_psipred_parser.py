"""Test functions for parsers.psipred_parser"""
import os
import unittest
from ample import constants
from ample.parsers import psipred_parser

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')

    def test_parse(self):
        ss2file = os.path.join(self.testfiles_dir, "1aba_.psipred_ss2")
        PA = psipred_parser.PsipredSs2Parser(ss2file)
        ref_ss2 = "CEEEEEECCCCCCCCHHHHHHHHHHHCCCCEEEEEECCCCCCCCHHHHHHHHHHHCCCCCCCCCCCEEEEECCEEEECHHHHHHHHC"
        ss2 = "".join([i.ss for i in PA.residues])
        self.assertEqual(ref_ss2, ss2 )
        
if __name__ == "__main__":
    unittest.main()
