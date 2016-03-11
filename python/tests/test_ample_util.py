"""Test functions for python.ample_util"""

import unittest
from ample.python import ample_util

class Test(unittest.TestCase):

    def test_ccp4Version(self):
        i, j, k = ample_util.ccp4_version()
        self.assertEqual((i, j, k),(7, 0, 4))
