"""Test functions for util.mrbump_util"""

import os
import unittest
from ample.util import mrbump_util

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-2 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')

    def test_process(self):
        with self.assertRaises(IOError):
            mrbdir = "/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_0/MRBUMP/MRBUMP"
            if not os.path.isdir(mrbdir): raise IOError
            rs = mrbump_util.ResultsSummary()
            print rs.summariseResults(mrbdir)
    
    def test_final_summary(self):
        with self.assertRaises(IOError):
            pkl = os.path.join(self.testfiles_dir, "resultsd.pkl")
            if not os.path.isfile(pkl): raise IOError
            with open(pkl) as f: d = cPickle.load(f)
            print mrbump_util.finalSummary(d)

if __name__ == "__main__":
    unittest.main()
