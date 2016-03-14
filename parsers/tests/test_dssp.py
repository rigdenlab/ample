
import os
import unittest
from ample.parsers import dssp_parser

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
        return

    def test_parse1(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        
        dsspLog = os.path.join(self.testfiles_dir,"2bhw.dssp")
        dsspP = dssp_parser.DsspParser( dsspLog )
        a = [' ', ' ', 'T', 'T', ' ', 'T', 'T', 'S', 'S', 'T', 'T', ' ', ' ', 
             ' ', 'T', 'T', 'G', 'G', 'G', ' ', ' ', 'S', ' ', ' ', 'T', 'T', 
             ' ', ' ', 'S', ' ', 'S', 'T', 'T', ' ', ' ', 'S', ' ', ' ', 'T', 
             'T', ' ', 'T', 'T', ' ', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'T', 
             'T', 'T', ' ', ' ', ' ', 'S', ' ', ' ', 'S', 'G', 'G', 'G', 'S', 
             'G', 'G', 'G', 'G', 'G', 'S', 'T', 'T', ' ', 'E', 'E', 'G', 'G', 
             'G', ' ', 'T', 'T', 'S', 'E', 'E', 'E', ' ', ' ', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'E', 'E', 'T', 'T', 'E', 'E', ' ', 'S', 
             'S', 'S', 'S', ' ', ' ', ' ', 'T', 'T', 'S', ' ', 'T', 'T', ' ', 
             'T', 'T', ' ', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', ' ', 'S', ' ', 
             'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', ' ', 'T', 'T', 
             'T', 'S', 'S', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', ' ', 'S', 
             ' ', ' ']
        self.assertEqual( dsspP.assignment[0], a)
        
    def test_parse2(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()

        dsspLog = os.path.join(self.testfiles_dir,"3ouf.dssp")
        dsspP = dssp_parser.DsspParser( dsspLog )

if __name__ == "__main__":
    unittest.main()
