"""Test functions for util.ample_util"""

import os
import sys
import unittest
from ample.util import ample_util

class Test(unittest.TestCase):

    def test_ccp4_version(self):
        i, j, k = ample_util.ccp4_version()
        self.assertEqual((i, j, 0),(7, 0, 0))
    
    def test_find_exe(self):
        gesamt_exe = os.path.basename(ample_util.find_exe("gesamt" + ample_util.EXE_EXT))
        self.assertTrue("gesamt.exe", gesamt_exe)
        spicker_exe = os.path.basename(ample_util.find_exe("spicker" + ample_util.EXE_EXT))
        self.assertTrue("spicker.exe", spicker_exe)
        theseus_exe = os.path.basename(ample_util.find_exe("theseus" + ample_util.EXE_EXT))
        self.assertTrue("theseus", theseus_exe)
        
if __name__ == "__main__":
    unittest.main()
