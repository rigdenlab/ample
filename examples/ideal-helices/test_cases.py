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
from test_funcs import AmpleException, parse_args

test_dict = {}

###############################################################################
#
# Ideal Helices
#
###############################################################################
args_ideal_helices =  [
                       '-fasta', '2OVC.fasta',
                       '-mtz', '2OVC-cad.mtz',
                       '-ideal_helices', 'True',
]

def test_ideal_helices(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not 'mrbump_results' in ad or not len(ad['mrbump_results']): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    if not ad['mrbump_results'][0]['SHELXE_CC'] > 25: raise AmpleException("SHELXE_CC criteria not met")
    return
        
test_dict['ideal_helices'] = { 'args' : args_ideal_helices,
                              'test' :  test_ideal_helices }

###############################################################################
#
# End Test Setup
#
###############################################################################

# Specify which directory these tests reside in
for name in test_dict.keys():
    test_dict[name]['directory'] = os.path.abspath(os.path.dirname(__file__))

if __name__ == '__main__':
    parse_args(test_dict)
