import os

if "CCP4" not in os.environ.keys():
    msg = "Cannot find CCP4 root directory"
    raise RuntimeError(msg)

__all__ = ["AMPLE_DIR", "SHARE_DIR", "AMPLE_CONFIG_FILE"]

AMPLE_DIR = os.path.join(os.environ["CCP4"], "lib", "python3.9", "site-packages", "ample")
SHARE_DIR = os.path.join(os.environ["CCP4"], "share", "ample")
AMPLE_CONFIG_FILE = os.path.join(SHARE_DIR, "include", "ample.ini")
AMPLE_LOGGER_CONFIG = os.path.join(SHARE_DIR, "include", "logging.json")
AMPLE_PKL = 'resultsd.pkl'
AMPLEDIR = 'AMPLE_'
I2DIR = 'AMPLEI2'
