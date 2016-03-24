#!/usr/bin/env ccp4-python

import argparse
import os

from ample.constants import SHARE_DIR
from ample.testing import test_funcs
from ample.testing.constants import EXAMPLE_DIRS, EXTRA_ARGS
from ample.testing.unittest_util import AMPLEUnittestFramework

def _integration(argd):
    all_test_cases = {}
    TEST_MODULE_NAME = 'test_cases'
    for directory in EXAMPLE_DIRS:
        test_module = test_funcs.load_module(TEST_MODULE_NAME, 
                                             [os.path.join(SHARE_DIR, "examples", 
                                                           directory)])
        for k, v in test_module.test_dict.iteritems():
            if k in all_test_cases:
                raise RuntimeError,"Duplicate key: {0}".format(k)
            all_test_cases[k] = v
    
    if argd['clean']:
        print "Cleaning test cases: {0}".format(all_test_cases.keys())
        test_funcs.clean(all_test_cases)
    else:
        print "Running test cases: {0}".format(all_test_cases.keys())
        test_funcs.run(all_test_cases, extra_args=EXTRA_ARGS, **argd)
    return

def _unittest(argd): 
    AMPLEUnittestFramework().run(buffer=argd['buffer'],
                                 cases=argd['test_cases'], 
                                 verbosity=argd['verbosity'])

def main():  
    desc = """ccp4-python -m ample.testing <command> [<args>]

Available tests include:
   integration     Integration testing of typical Ample routines
   unittest        Unittesting of all Ample subroutines
"""
    
    parser = argparse.ArgumentParser(prog="run_tests.py",
                                     usage=desc)
    suboptions = parser.add_subparsers(help="Testing framework options")
    
    # Integration testing using examples
    integ = suboptions.add_parser("integration", 
                                  help="Integration testing with examples")
    integ.set_defaults(which="integration")
    integ.add_argument('-c', '--clean', action='store_true', default=False,
                       help="Clean up all test files/directories")
    integ.add_argument('-n', '--nproc', type=int, default=1,
                       help="Number of processors to run on (1 per job)")
    integ.add_argument('-d', '--dry_run', action='store_true', default=False,
                       help="Don\'t actually run the jobs")
    integ.add_argument('-r', '--rosetta_dir',
                       help="Location of rosetta installation directory")
    integ.add_argument('-s', '--submit_cluster', action='store_true', default=False,
                       help="Submit to a cluster queueing system")
    integ.add_argument('test_cases', nargs='*',
                       help="A list of test cases to run")
    
    # Function unittesting
    unit = suboptions.add_parser("unittest", 
                                 help="Unittest all functions")
    unit.set_defaults(which="unittest")
    unit.add_argument('-b', dest='buffer', action="store_false", default=True,
                      help="debugging by printing print messages")
    unit.add_argument('test_cases', nargs='*',
                      help="A list of test cases to run")
    unit.add_argument('-v', dest="verbosity", default=2, type=int, 
                      help="level of verbosity [default: 2]")
    
    argd = vars(parser.parse_args())
      
    if argd['which'] == "integration" :
        _integration(argd)
    elif argd['which'] == 'unittest':
        _unittest(argd)
    
if __name__ == "__main__":
    main()
