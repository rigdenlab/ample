from enum import Enum
import json
import logging.config
import os
import sys

from ample.constants import AMPLE_LOGGER_CONFIG


class LogColors(Enum):
    """Color container for log messages"""
    CRITICAL = 31
    DEBUG = 34
    DEFAULT = 0
    ERROR = 31
    WARNING = 33
 
 
class LogColorFormatter(logging.Formatter):
    """Formatter for log messages"""
    def format(self, record):   
        if record.levelname in LogColors.__members__:
            prefix = '\033[1;{}m'.format(LogColors[record.levelname].value)
            postfix = '\033[{}m'.format(LogColors["DEFAULT"].value)
            record.msg = os.linesep.join([prefix + msg + postfix for msg in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)
 
 
def setup_logging(argso):
    """Read JSON config for logger and return root logger
    
    Also sets the path to the AMPLE logfile in the dictionary (required for pyrvapi)"""
    if not os.path.isfile(AMPLE_LOGGER_CONFIG):
        raise RuntimeError("Cannot find AMPLE_LOGGER_CONFIG file: {}".format(AMPLE_LOGGER_CONFIG))
    with open(AMPLE_LOGGER_CONFIG, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
    try:
        argso['ample_log'] = os.path.abspath(config['handlers']['file_handler']['filename'])
    except KeyError:
        argso['ample_log'] = None
    return logging.getLogger()   


def setup_console_logging(level=logging.INFO,
                          formatstr='%(message)s\n'):
    """
    Set up logging to the console - required for the individual modules.
    
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
    Set up logging to a file - required for the individual modules.
    
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