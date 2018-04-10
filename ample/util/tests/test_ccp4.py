"""Test functions for util.ccp4"""

import unittest
from ample.util import ccp4

class Test(unittest.TestCase):

    def test_version_1(self):
        self.assertTrue(ccp4.CCP4Version().version > (7, 0))


if __name__ == "__main__":
    unittest.main()
