#!/usr/bin/env ccp4-python
"""
This is AMPLE

This script is named ample.py due to a problem with running the multiprocessing (which is used to parallelise
the running of jobs on a local machine - see python/workers.py) module under windows.
The multiprocessing module on windows requires that it can import the main module, and the import machinery
requires that any file being imported is named <foo>.py, and any changes to this would require hacking the 
multiprocessing module, so to avoid this, our script must be called ample.py
"""
import os
import sys

# Add the ample python folder to the PYTHONPATH
if "CCP4_AMPLE_ROOT" in os.environ.keys() and "CCP4" in os.environ.keys():
    root = os.environ["CCP4_AMPLE_ROOT"]
elif "CCP4" in os.environ.keys():
    root = os.path.join(os.environ["CCP4"], "share", "ample")
else:
    raise RuntimeError('CCP4 not found')

sys.path.append( os.path.join( root, "python" ) )

# our imports
import ample_main # Order is important as ample_main sets up the root console logger
import ample_exit

# Run ample
try:
    ample_main.Ample().run()
except Exception as e:
    msg = "Error running main AMPLE program: {0}".format(e.message)
    ample_exit.exit_error(msg, sys.exc_info()[2])

