import logging
import sys

def setup_console_logging(level=logging.INFO,
                          formatstr='%(message)s\n' # Always add a blank line after every print
                          ):
    """
    Set up logging to the console for the root logger.
    
    Parameters
    ----------
    level : int
        Sets the threshold for the console output to level.
    formatstr : str
        The string used to format the log messages
        
    Returns
    -------
    logger : :obj:logging.logger
        The root logger
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # Seems they changed the api in python 2.6->2.7
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(stream=sys.stdout)
    cl.setLevel(level)
    formatter = logging.Formatter(formatstr) 
    cl.setFormatter(formatter)
    logger.addHandler(cl)
    
    return logger
    
def setup_file_logging(logfile,
                       level=logging.DEBUG,
                       formatstr='%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
    """
    Set up logging to a file for the root logger.
    
    Parameters
    ----------
    logfile : str
        The path to the logfile that output will be written to.
    level : int
        Sets the threshold for the console output to level.
    formatstr : str
        The string used to format the log messages
        
    Returns
    -------
    logger : :obj:logging.logger
        The root logger
    """
    logger = logging.getLogger()
    fl = logging.FileHandler(logfile)
    fl.setLevel(level)
    formatter = logging.Formatter(formatstr)
    fl.setFormatter(formatter)
    logger.addHandler(fl)
    
    return logger
