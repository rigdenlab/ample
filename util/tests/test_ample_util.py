"""Test functions for util.ample_util"""

import os
import sys
import unittest
from ample.util import ample_util

class Test(unittest.TestCase):

    def test_ccp4_version(self):
        i, j, k = ample_util.ccp4_version()
        self.assertEqual((i, j, 0),(7, 0, 0))

    @unittest.skipIf(sys.platform.startswith("win"), "")
    def test_find_exe_UNIX(self):
        self.assertTrue("gesamt", os.path.basename(ample_util.find_exe("gesamt")))
        self.assertTrue("spicker", os.path.basename(ample_util.find_exe("spicker")))
        self.assertTrue("theseus", os.path.basename(ample_util.find_exe("theseus")))
        
    @unittest.skipUnless(sys.platform.startswith("win"), "")
    def test_find_exe_WINDOWS(self):
        self.assertTrue("gesamt.exe", os.path.basename(ample_util.find_exe("gesamt")))
        self.assertTrue("spicker.exe", os.path.basename(ample_util.find_exe("spicker")))
        self.assertTrue("theseus.exe", os.path.basename(ample_util.find_exe("theseus")))
        
if __name__ == "__main__":
    unittest.main()