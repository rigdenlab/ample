#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''
import os
import sys

sys.path.append('/opt/ample-dev1/python')
import test_funcs


dirs = [
    'homologs',
    'missing-domain.1k04',
    'toxd-example',
    'transmembrane.3LBW'
    ]

dirs = [
    'ideal-helices',
    'nmr.remodel',
    'nmr.truncate',
    'toxd-example'
]

EXTRA_ARGS = [ '-no_gui','True',
              # '-do_mr','False',
               ]

# Get the command-line arguments
argd = test_funcs.parse_args()

owd = os.getcwd()
all_test_cases = {}
for d in dirs:
    d = os.path.abspath(d)
    # possibly a bit clunky - we add the test directory to the path so we can import the
    # test dict and them remove it from the sys.path so we get the next module next time
    sys.path.append(d)
    import test_cases
    for k, v in test_cases.test_dict.iteritems():
        if k in all_test_cases:
            raise RuntimeError,"Duplicate key: {0}".format(k)
        all_test_cases[k] = v
        
    # clean path and unload module
    sys.path.remove(d)
    del sys.modules['test_cases']


if argd['clean']:
    print "Cleaning test cases: {0}".format(all_test_cases.keys())
    test_funcs.clean(all_test_cases)
else:
    print "Running test cases: {0}".format(all_test_cases.keys())
    test_funcs.run(all_test_cases, extra_args=EXTRA_ARGS, **argd)


