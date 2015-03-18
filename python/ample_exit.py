'''
Created on Mar 18, 2015

@author: jmht
'''

import logging
import sys
import traceback

# external imports
try: import pyrvapi
except: pyrvapi=None

def exit(msg):
    logger = logging.getLogger()
    
    #header="**** AMPLE ERROR ****\n\n"
    header="*"*70+"\n"
    header+="*"*20 + " "*10 + "AMPLE ERROR" + " "*10 +"*"*19 + "\n" 
    header+="*"*70+"\n\n"
    
    footer="\n\n" + "*"*70+"\n\n"
    
    # Bit dirty - get the name of the debug log file
    debug_log=None
    for d in logger.handlers:
        n='baseFilename'
        if hasattr(d,n) and d.level==logging.DEBUG:
            debug_log=getattr(d, n)
    if debug_log:
        footer+="More information may be found in the debug log file: {0}\n".format(debug_log) 
    
    footer += "\nIf you believe that this is an error with AMPLE, please email: ccp4@stfc.ac.uk\n"
    footer += "providing as much information as you can about how you ran the program.\n"
    if debug_log: 
        footer += "\nPlease include the debug logfile with your email: {0}\n".format(debug_log)   
    
    # String it all together
    msg=header + msg + footer
    
    logger.critical(msg)
    
    # Get traceback of where we failed for the log file
    logger.debug("AMPLE EXITING AT...")
    logger.debug("".join(traceback.format_list(traceback.extract_stack())))
    
    # Make sure the error widget is updated
    if pyrvapi: pyrvapi.rvapi_flush()
    
    sys.exit(1)
  
                 
