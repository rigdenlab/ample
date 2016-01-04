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

EXTRA_ARGS = ['-no_gui','True']

argd = test_funcs.parse_args()

owd = os.getcwd()
for d in dirs:
    d = os.path.abspath(d)
    # possibly a bit clunky - we add the test directory to the path so we can import the
    # test dict and them remove it from the sys.path so we get the next module next time
    sys.path.append(d)
    import test_cases
    os.chdir(d)
    if argd['clean']:
        print "Cleaning directory: {0}".format(d)
        test_funcs.clean(test_cases.test_dict)
    else:
        print "RUNNING TESTS IN DIRECTORY: {0}".format(d)
        test_funcs.run(test_cases.test_dict, extra_args=EXTRA_ARGS, **argd)

    # clean path and unload module
    os.chdir(owd)
    sys.path.remove(d)
    del sys.modules['test_cases']

