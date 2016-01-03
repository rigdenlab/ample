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
    'toxd-example',
]

EXTRA_ARGS = ['-no_gui','True']

argd = test_funcs.parse_args()

owd = os.getcwd()
for d in dirs:
    # possibly a bit clunky - we add the test directory to the path so we can import the
    # test dict and them remove it from the sys.path
    sys.path.append(d)
    from test_cases import test_dict
    os.chdir(d)
    if argd['clean']:
        test_funcs.clean(test_dict)
    else:
        print "RUNNING TESTS IN DIRECTORY: {0}".format(os.path.abspath(d))
        test_funcs.run(test_dict, extra_args=EXTRA_ARGS, **argd)
    os.chdir(owd)
    sys.path.pop()
    

