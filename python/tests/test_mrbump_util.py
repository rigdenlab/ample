"""Test functions for python.ample_mrbump"""

import unittest
from ample.python import mrbump_util

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
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)
        return

    def XtestProcess(self):
        """Parse a results file"""
        
        mrbdir = "/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_0/MRBUMP/MRBUMP"
        rs = mrbump_util.ResultsSummary()
        print rs.summariseResults(mrbdir)
        return
    
    def test_finalSummary(self):
        """Parse a results file"""
        pkl = os.path.join(self.testfiles_dir, "resultsd.pkl")
        with open(pkl) as f: d = cPickle.load(f)
        print mrbump_util.finalSummary(d)
        return