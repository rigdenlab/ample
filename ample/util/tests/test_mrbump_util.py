"""Test functions for util.mrbump_util"""

import cPickle
import os
import unittest

from ample import constants
from ample.util import mrbump_util

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')

    def test_process(self):
        mrbdir = os.path.join(self.testfiles_dir, "MRBUMP")
        if not os.path.isdir(mrbdir): return
        rs = mrbump_util.ResultsSummary()
        print rs.summariseResults(mrbdir)

    def test_final_summary(self):
        pkl = os.path.join(self.testfiles_dir, "resultsd.pkl")
        if not os.path.isfile(pkl): return
        with open(pkl) as f: d = cPickle.load(f)
        print mrbump_util.finalSummary(d)

if __name__ == "__main__":
    unittest.main()
