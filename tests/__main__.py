"""main routine for Ample testing"""
import sys
from ample.python import ample_exit
from ample.tests import run_tests

try:
    run_tests.main()
except Exception as e:
    msg = "Error running Ample testsuite: {0}".format(e.message)
    ample_exit.exit_error(msg, sys.exc_info()[2])
