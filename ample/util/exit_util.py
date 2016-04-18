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

def exit_error(msg, ample_tb=None):
    """Exit on error collecting as much information as we can.
    
    args:
    ample_tb - this can be got from sys.exc_info()[2]
    """
    # Get the root logger 
    logger = logging.getLogger()
    
    # An error may have occured before we started logging so we need to create one here
    if not logger.handlers:
        logging.basicConfig(format='%(message)s\n', level=logging.DEBUG)
        logger = logging.getLogger()
    
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
    
    # Get traceback of where we failed for the log file
    if not ample_tb:
        ample_tb = traceback.extract_stack()
    else:
        ample_tb = traceback.extract_tb(ample_tb)
    
    msg = "AMPLE EXITING AT..." + os.linesep + "".join(traceback.format_list(ample_tb))
    if debug_log:
        logger.debug(msg)
    else:
        # If we don't have a debug file we want to output the traceback to the console
        logger.info(msg)
    
    # Make sure the error widget is updated
    if pyrvapi: pyrvapi.rvapi_flush()
    
    sys.exit(1)
  
                 
