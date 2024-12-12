#!/usr/bin/env ccp4-python

import argparse

from ample.testing import integration_util, unittest_util
from ample.util.argparse_util import add_core_options, add_cluster_submit_options

__author__ = "Felix Simkovic"
__date__ = "25-Mar-2016"


def run_integration(argd):
    m = integration_util.AMPLEIntegrationFramework(test_cases=argd['test_cases'], run_dir=argd['run_dir'])
    if argd['clean']:
        m.clean() if argd['test_cases'] else m.clean(clean_dir=True)
    else:
        m.run(**argd)


def run_unittest(argd):
    m = unittest_util.AMPLEUnittestFramework()
    m.run(buffer=argd['buffer'], cases=argd['test_cases'], verbosity=argd['verbosity'])


def main():
    desc = """ccp4-python -m ample.testing <command> [<args>]

Available tests include:
   integration     Integration testing of typical Ample routines
   unittest        Unittesting of all Ample subroutines
"""

    parser = argparse.ArgumentParser(prog="run_tests.py", usage=desc)

    # Add options that work on all runtypes
    add_core_options(parser)

    suboptions = parser.add_subparsers(help="Testing framework options")

    integ = suboptions.add_parser("integration", help="Integration testing with examples")
    integ.set_defaults(which="integration")
    integration_util.add_cmd_options(integ)
    add_cluster_submit_options(integ)

    unit = suboptions.add_parser("unittest", help="Unittest all functions")
    unit.set_defaults(which="unittest")
    unittest_util.add_cmd_options(unit)

    argd = vars(parser.parse_args())

    which_test = argd['which']
    if which_test == 'integration':
        run_integration(argd)
    elif which_test == 'unittest':
        run_unittest(argd)


if __name__ == "__main__":
    main()
