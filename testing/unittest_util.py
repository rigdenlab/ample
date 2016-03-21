# Module containing a framework for unittesting of AMPLE modules. This framework 
# serves the purpose of loading all/selected test cases and running them.
#
# Author: hlfsimko
# Date: 21-Mar-2016

import os
from ample.testing.constants import AMPLE_DIR
from unittest import TestLoader, TextTestRunner, TestSuite

class AMPLEUnittestFramework(object):
    """Framework to run Ample unittesting"""
    
    def run(self, buffer=True, cases=None, pattern="test*.py", verbosity=2):
        """main routine for running the test cases
        
        :buffer: display STDOUT messages
        :case: list of folders to be tested
        :pattern: test file name pattern
        :verbosity: chatty output
        """
        if cases:
            suite = self._load_for_subselection(cases, pattern)
        else:
            suite = TestLoader().discover(AMPLE_DIR, pattern=pattern,
                                          top_level_dir=AMPLE_DIR)
        TextTestRunner(verbosity=verbosity, buffer=buffer).run(suite)

    def _load_for_subselection(self, cases, pattern):
        suite = TestSuite()
        for case in cases:
            path = os.path.join(AMPLE_DIR, case)
            _suite = TestLoader().discover(path, pattern=pattern,
                                           top_level_dir=AMPLE_DIR)
            suite.addTests(_suite)
            del _suite
        return suite