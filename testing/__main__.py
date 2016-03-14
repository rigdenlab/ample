"""main routine for Ample testing"""
import sys
from ample.util import exit_util
from ample.testing import run_tests

try:
    run_tests.main()
except Exception as e:
    msg = "Error running Ample testsuite: {0}".format(e.message)
    exit_util.exit_error(msg, sys.exc_info()[2])
