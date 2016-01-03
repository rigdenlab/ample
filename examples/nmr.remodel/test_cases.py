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
# NMR Remodelling
#
###############################################################################
args =  [
    '-name', '1t00',
    '-fasta', '1T00.fasta',
    '-mtz', '1t00.mtz',
    '-rosetta_dir', '/opt/rosetta-3.5',
    '-nmr_model_in', '2DIZ.pdb',
    '-nmr_remodel', 'True',
    '-frags_3mers', '1t00.200.3mers',
    '-frags_9mers', '1t00.200.9mers',
    '-nmr_process', '1',
]

def testf(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ('mrbump_results' in ad or len(ad['mrbump_results'])): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    if not ad['mrbump_results'][0]['SHELXE_CC'] > 25: raise AmpleException("SHELXE_CC criteria not met")
    return
        
test_dict['nmr_remodel'] = { 'args' : args,
                              'test' :  testf }

###############################################################################
#
# End Test Setup
#
###############################################################################

if __name__ == '__main__':
    test_funcs.parse_args(test_dict)
