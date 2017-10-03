"""Test functions for util.ample_util"""

import cPickle
import os
import unittest
from ample.util import ample_util
from ample import constants

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
        
    def test_resultsd_fix_path(self):
        with open(os.path.join(constants.SHARE_DIR,'testfiles','resultsd.pkl')) as f: optd = cPickle.load(f)
        d = ample_util.resultsd_fix_path(optd, newroot='/foo/bar')
        self.assertEqual(d['fasta'],'/foo/bar/AMPLE_0/ampl_.fasta')
        self.assertEqual(d['mrbump_results'][0]['PHASER_logfile'],'/foo/bar/AMPLE_0/MRBUMP/search_c1_t69_r3_polyAla_mrbump/data/loc0_ALL_c1_t69_r3_polyAla/unmod/mr/phaser/phaser_loc0_ALL_c1_t69_r3_polyAla_UNMOD.log')
        
if __name__ == "__main__":
    unittest.main()
