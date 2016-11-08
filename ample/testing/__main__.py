"""main routine for Ample testing"""
import logging
import os
import sys
from ample.util import exit_util
from ample.testing import run_tests

logging.basicConfig()
logger = logging.getLogger()
#jmht - the default handler is set to DEBUG so we need to set it to CRITICAL
# to avoid drowning in output
logger.handlers[0].setLevel(logging.CRITICAL)

#############################################################################
## Multiprocessing crashes on Windows when running multiple jobs.
#  Issue recorded 
#       1) https://docs.python.org/2/library/multiprocessing.html#windows
if sys.platform.startswith("win"):
    msg = """
*****************************************************************************
A bug prevents you from invoking our testing framework via the module loader. 
                                                                              
Please invoke using the following command:                                    
                                                                              
% ccp4-python {0}{1}run_tests.py <command> [<args>]                           
*****************************************************************************
"""
    msg.format(os.path.dirname(__file__), os.sep)
    logger.critical(msg)
    sys.exit(1)

#############################################################################
## On Unix systems we can run as normal
try:
    run_tests.main()
except Exception as e:
    msg = "Error running Ample testsuite: {0}".format(e.message)
    logger.critical(msg)
    exit_util.exit_error(msg, sys.exc_info()[2])
