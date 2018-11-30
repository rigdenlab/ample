
import argparse
import sys

if sys.version_info.major < 3:
    from urllib2 import urlopen, HTTPError
else:
    from urllib.error import URLError
    from urllib.request import urlopen

from ample.testing import integration_util
from ample.testing import run_tests
from ample.util import ample_util

__author__ = "Felix Simkovic and Jens Thomas"
__date__ = "26-Mar-2016"

def found_exe(e):
    try:
        ample_util.find_exe(e)
    except:
        return False
    return True

def internet_on():
    # Taken from stackoverflow.com/questions/3764291/checking-network-connection
    try:
        reponse = urlopen("http://google.com", timeout=5)
    except URLError:
        return False
    return True

def parse_args(test_dict):
    """wrapper function for running test cases from scripts directly"""
    parser = argparse.ArgumentParser()
    integration_util.add_cmd_options(parser)
    argd =  vars(parser.parse_args())
    if argd['test_cases']:
        argd['test_cases'] = list(set(test_dict.keys()) & set(argd['test_cases']))
    else:
        argd['test_cases'] = test_dict.keys()
    run_tests.run_integration(argd)
    
    
