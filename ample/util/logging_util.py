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
    with open(AMPLE_LOGGER_CONFIG, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
    argso['ample_log'] = os.path.abspath(config['handlers']['file_handler']['filename'])
    return logging.getLogger()   


def setup_console_logging(level=logging.INFO,
                          formatstr='%(message)s\n'):
    """
    Set up logging to the console - required for the test framework.
    
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