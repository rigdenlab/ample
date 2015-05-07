#!/usr/bin/env ccp4-python
"""Script to run all the AMPLE tests"""
import os
import sys
import unittest

###############################################################
#
# Add the python directory to the path
#
###############################################################
thisd =  os.path.abspath(os.path.dirname(__file__))
paths = thisd.split(os.sep)
ample_dir = os.sep.join(paths[:-1])
python_dir = os.path.join(ample_dir,'python')
sys.path.insert(0,python_dir)

loader = unittest.defaultTestLoader

tests = loader.discover(python_dir, pattern="*.py")

unittest.TextTestRunner(verbosity=2).run(tests) 
