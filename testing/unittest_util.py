
# Module containing a framework for unittesting of AMPLE modules. 
# This framework serves the purpose of loading and running all/selected test cases.
#
# Author: hlfsimko
# Date: 21-Mar-2016

import glob
import os
from ample.testing.constants import AMPLE_DIR
from unittest import TestLoader, TextTestRunner, TestSuite

class AMPLEUnittestFramework(object):
    """Framework to run Ample unittesting"""
    
    def run(self, buffer=False, cases=None, pattern="test*.py", verbosity=2):
        """main routine for running the test cases"""
        
        suite = SuiteLoader().load_suite(AMPLE_DIR, cases=cases, pattern=pattern)
        
        if int(suite.countTestCases()) > 0:
            TextTestRunner(verbosity=verbosity, buffer=buffer).run(suite)
            
        return

class SuiteLoader(object):
    """Loader designed to obtain all test cases in a package"""
    
    def load_suite(self, directory, pattern="test*.py", cases=None):
        """function to load a unittest test suite"""
        # If we do not have any test cases then we can search for some in
        # the specified directory.
        if not cases:
            search_pattern = os.path.join(directory, "*")
            cases = [ os.path.basename(folder) for folder in \
                        glob.glob(search_pattern) if os.path.isdir(folder) ]
            
        suite = self._load_suite(cases, pattern, directory)
        return suite

    def _load_suite(self, cases, pattern, directory):
        suite = TestSuite()
        for case in cases:
            path = os.path.join(directory, case)
            
            try:
                _suite = TestLoader().discover(path, pattern=pattern,
                                               top_level_dir=directory)
                suite.addTests(_suite)
                del _suite
            except ImportError:
                print "*** no test cases found for: {0} ***".format(case)
            
        return suite