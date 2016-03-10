"""Facility to run all unittests"""

import unittest

if __name__ == "__main__":
    suite = unittest.TestLoader().discover("..", pattern="test*.py")
    unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)
