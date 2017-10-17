'''
Created on Mar 18, 2015

@author: jmht
'''

import logging
import os
import sys
import traceback

# external imports
try: import pyrvapi
except: pyrvapi=None

def _debug_logfile(logger):
    debug_log = None
    if logger.handlers:
        for d in logger.handlers:
            n='baseFilename'
            if hasattr(d,n) and d.level==logging.DEBUG:
                debug_log=getattr(d, n)
    return debug_log

def exit_error(*args, **kwargs):
    """Exit on error collecting as much information as we can.
    
    Parameters
    ----------
    message : str, optional
      A error message to print
    
    Notes
    -----
    This previously accepted two arguments of a string to print as an error message and
    an exception traceback.
    We now just use sys.exch_info() so the messsage argument is no longer required but optional.
    While we refactor the code, we'll use *args to get any argument parameters passed in
    as the first argument.

    """
    # Get the root logger 
    logger = logging.getLogger()
    
    # An error may have occured before we started logging so we need to create one here
    if not logger.handlers:
        logging.basicConfig(format='%(message)s\n', level=logging.DEBUG)
        logger = logging.getLogger()

    
    exc_type, exc_value, exc_traceback = sys.exc_info()
    msg = kwargs.get('message')
    if msg is None:
        if len(args) >= 1: # Fix for old cases
            msg = args[0]
        else:
            msg = "{0}: {1}".format(exc_type.__name__, exc_value.message)

    #header="**** AMPLE ERROR ****\n\n"
    header="*"*70 + "\n"
    header+="*"*20 + " "*10 + "AMPLE ERROR" + " "*10 +"*"*19 + "\n" 
    header+="*"*70 + "\n\n"
    
    # Create the Footer 
    footer="\n\n" + "*"*70+"\n\n"
    
    # Get the name of the debug log file
    debug_log = _debug_logfile(logger)
    if debug_log: footer += "More information may be found in the debug log file: {0}\n".format(debug_log) 
    
    footer += "\nIf you believe that this is an error with AMPLE, please email: ccp4@stfc.ac.uk\n"
    footer += "providing as much information as you can about how you ran the program.\n"
    if debug_log: footer += "\nPlease include the debug logfile with your email: {0}\n".format(debug_log)   
    
    # String it all together
    msg = header + msg + footer
    
    # Print out main message
    logger.critical(msg)
    
    # If we were called without an exception being raised, we just print the current stack
    if exc_traceback is None: exc_traceback = traceback.extract_stack()
    msg = "AMPLE EXITING AT..." + os.linesep + "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    if debug_log:
        logger.debug(msg)
    else:
        # If we don't have a debug file we want to output the traceback to the console
        logger.info(msg)
    
    # Make sure the error widget is updated
    if pyrvapi: pyrvapi.rvapi_flush()
    
    sys.exit(1)
  
                 
