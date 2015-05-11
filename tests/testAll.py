#!/usr/bin/env ccp4-python
"""Script to run all the AMPLE tests"""
import glob
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
# Don't use discover as we don't want to run all tests
#tests = loader.discover(python_dir, pattern="*.py")
skip = ['benchmark', 'clusterize', 'mrbump_results', 'rosetta_model']
test_suite = unittest.TestSuite()
# Get list of modules that are not in skip and add the test sets to the test suite
for m in glob.glob(os.path.join(python_dir,"*.py")):
    m = os.path.splitext(os.path.basename(m))[0]
    if m not in skip:
        m = __import__(m)
        test_suite.addTests(loader.loadTestsFromModule(m))

unittest.TextTestRunner(verbosity=2).run(test_suite) 
