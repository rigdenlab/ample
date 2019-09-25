"""main routine for Ample testing"""
import logging
import os
import sys
from ample.util import exit_util, logging_util
from ample.testing import run_tests

logger = logging_util.setup_console_logging(
    level=logging.INFO, formatstr='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

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
