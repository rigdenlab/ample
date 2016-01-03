#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import cPickle
import os
import sys

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
import test_funcs

test_dict = {}

###############################################################################
#
# NMR Truncation
#
###############################################################################
args =  [
    '-name', '102l',
    '-fasta', '102l.fasta',
    '-mtz', '102L.mtz',
    '-nmr_model_in', '2LC9.pdb',
]

def testf(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ('mrbump_results' in ad or len(ad['mrbump_results'])): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    if not ad['mrbump_results'][0]['SHELXE_CC'] > 25: raise AmpleException("SHELXE_CC criteria not met")
    return
        
test_dict['nmr_truncate'] = { 'args' : args,
                              'test' :  testf }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(test_dict)
