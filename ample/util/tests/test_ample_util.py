"""Test functions for util.ample_util"""

import pickle
import os
import sys
import unittest
from ample.util import ample_util
from ample.constants import AMPLE_PKL, SHARE_DIR


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        cls.testfiles_dir = os.path.join(SHARE_DIR, 'testfiles')

    def test_find_exe(self):
        gesamt_exe = os.path.basename(ample_util.find_exe("gesamt" + ample_util.EXE_EXT))
        self.assertEqual("gesamt" + ample_util.EXE_EXT, gesamt_exe)
        spicker_exe = os.path.basename(ample_util.find_exe("spicker" + ample_util.EXE_EXT))
        self.assertEqual("spicker" + ample_util.EXE_EXT, spicker_exe)
        theseus_exe = os.path.basename(ample_util.find_exe("theseus" + ample_util.EXE_EXT))
        self.assertEqual("theseus" + ample_util.EXE_EXT, theseus_exe)

    def test_resultsd_fix_path(self):
        with open(os.path.join(SHARE_DIR, 'testfiles', AMPLE_PKL), 'rb') as f:
            if sys.version_info.major == 3:
                optd = pickle.load(f, encoding='latin1')
            else:
                optd = pickle.load(f)
        d = ample_util.amoptd_fix_path(optd, newroot='/foo/bar')
        self.assertEqual(d['fasta'], '/foo/bar/ampl_.fasta')
        self.assertEqual(
            d['mrbump_results'][0]['PHASER_logfile'],
            '/foo/bar/MRBUMP/search_c1_t69_r3_polyAla_mrbump/data/loc0_ALL_c1_t69_r3_polyAla/unmod/mr/phaser/phaser_loc0_ALL_c1_t69_r3_polyAla_UNMOD.log',
        )

    def test_extract_tar(self):
        tarfile = os.path.join(self.testfiles_dir, 'result.tar.bz2')

        # Extract single file
        filename = 'alldecoy.pdb'
        files = ample_util.extract_tar(tarfile, filenames=[filename])
        self.assertEqual(len(files), 1)
        fpath = os.path.abspath(filename)
        self.assertEqual(fpath, files[0])
        os.unlink(fpath)

        # Extract multiple files
        suffixes = ['.pdb']
        files = ample_util.extract_tar(tarfile, suffixes=suffixes)
        self.assertEqual(len(files), 6)
        for f in files:
            os.unlink(f)


if __name__ == "__main__":
    unittest.main()
