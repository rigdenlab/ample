"""Worker functions for job execution"""

from __future__ import print_function

__author__ = "Jens Thomas"
__date__ = "28 Feb 2013"
__version__ = "1.0"

import multiprocessing
import os
import sys

from ample.util import ample_util


def worker(inqueue, early_terminate=False, check_success=None):
    """Worker process to run MrBump jobs until no more left.

    This function keeps looping over the inqueue, removing jobs from the 
    inqueue until there are no more left. It checks if a jobs has succeeded
    and if so it will terminate.

    Parameters
    ----------
    inqueue : :obj:`Queue`
       A Python Queue object
    early_terminate : bool
       Terminate on first success or continue running
    check_success : callable
       A callable to check the success status of a job
    
    Warnings
    --------
    This needs to import the main module that it lives in so maybe this should
    live in a separate module?

    """
    if early_terminate:
        assert callable(check_success)

    success = True
    while True:
        if inqueue.empty():
            print ("worker {0} got empty inqueue".format(multiprocessing.current_process().name))
            rcode = 0 if success else 1
            sys.exit(rcode)

        # Got a script so run
        job = inqueue.get()

        # Get name from script
        print ("Worker {0} running job {1}".format(multiprocessing.current_process().name, job))
        directory, sname = os.path.split(job)
        jobname = os.path.splitext(sname)[0]

        # Change directory to the script directory
        os.chdir(directory)
        retcode = ample_util.run_command([job], logfile=jobname + ".log", dolog=False, check=True)

        # Can we use the retcode to check?
        # REM - is retcode object
        if retcode != 0:
            print ("WARNING! Worker {0} got retcode {1}".format(multiprocessing.current_process().name, retcode))
            success = False

        # Now check the result if early terminate
        if early_terminate:
            if check_success(job):
                print ("Worker {0} job succeeded".format(multiprocessing.current_process().name))
                sys.exit(0)

    print ("worker {0} FAILED!".format(multiprocessing.current_process().name))
    sys.exit(1)
