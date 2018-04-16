"""Test functions for util.tm_util"""

import unittest
from ample.testing import test_funcs
from ample.util import ample_util, tm_util

@unittest.skipUnless(test_funcs.found_exe("TMscore" + ample_util.EXE_EXT), "TMscore exec missing")
class TestTM(unittest.TestCase):

    def test_gaps_1(self):
        gaps = tm_util.TMscore("TMscore", wdir=".")._find_gaps("AAAA---AA--AA")
        ref_gaps = [False, False, False, False, True, True, True, False, False, True, True, False, False]
        self.assertEqual(ref_gaps, gaps)

    def test_gaps_2(self):
        gaps = tm_util.TMscore("TMscore", wdir=".")._find_gaps("---AA--AA")
        ref_gaps = [True, True, True, False, False, True, True, False, False]
        self.assertEqual(ref_gaps, gaps)

    def test_gaps_3(self):
        gaps = tm_util.TMscore("TMscore", wdir=".")._find_gaps("-AAA--")
        ref_gaps = [True, False, False, False, True, True]
        self.assertEqual(ref_gaps, gaps)

if __name__ == "__main__":
    unittest.main()
