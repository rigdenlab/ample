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

    def test_final_summary(self):
        pkl = os.path.join(self.testfiles_dir, "resultsd.pkl")
        if not os.path.isfile(pkl): return
        with open(pkl) as f: d = cPickle.load(f)
        summary = mrbump_util.finalSummary(d)
        self.assertIsNotNone(summary)
        
    def test_topfiles(self):
        topf = mrbump_util.ResultsSummary(results_pkl=os.path.join(self.testfiles_dir, "resultsd.pkl")).topFiles()
        self.assertEqual(len(topf),3)
        self.assertIn('info',topf[2])

if __name__ == "__main__":
    unittest.main()
