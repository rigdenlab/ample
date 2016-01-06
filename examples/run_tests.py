#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''
import os
import sys

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-1 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

dirs = [
    'homologs',
    'missing-domain.1k04',
    'toxd-example',
    'transmembrane.3LBW'
    ]

# List of which test directories to process
dirs = [
    'ideal-helices',
    'nmr.remodel',
    'nmr.truncate',
    'toxd-example'
]

# Any args that are to be added/updated
# They *MUST* be paired - i.e. nothing with more then one option
EXTRA_ARGS = [ '-no_gui','True',
              # '-do_mr','False',
               ]

# Get the command-line arguments
argd = test_funcs.parse_args()

all_test_cases = {}
TEST_MODULE_NAME = 'test_cases'
for directory in dirs:
    test_module = test_funcs.load_module(TEST_MODULE_NAME, [os.path.abspath(directory)])
    for k, v in test_module.test_dict.iteritems():
        if k in all_test_cases:
            raise RuntimeError,"Duplicate key: {0}".format(k)
        all_test_cases[k] = v

if argd['clean']:
    print "Cleaning test cases: {0}".format(all_test_cases.keys())
    test_funcs.clean(all_test_cases)
else:
    print "Running test cases: {0}".format(all_test_cases.keys())
    test_funcs.run(all_test_cases, extra_args=EXTRA_ARGS, **argd)


