#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''

import cPickle
import glob
import os
import sys

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-2 ])
sys.path.append(os.path.join(AMPLE_DIR,'python'))
from test_funcs import AmpleException, parse_args

test_dict = {}

# vanilla test
args_vanilla = [
    '-fasta', 'toxd_.fasta',
    '-mtz', '1dtx.mtz',
    '-percent', '50',
]

###############################################################################
#
# Making models
#
###############################################################################

args_rosetta_modelling = args_vanilla + [
    '-rosetta_dir', '/opt/rosetta-3.5',
    '-frags_3mers', 'aat000_03_05.200_v1_3',
    '-frags_9mers', 'aat000_09_05.200_v1_3',
    '-nmodels', '30',
]

def test_rosetta_modelling(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ad['models_dir'] or not len(glob.glob(ad['models_dir']+"/*.pdb")) == 30: raise AmpleException("Incorrect number of models")
    if not ('ensembles' in ad or len(ad['ensembles'])): raise AmpleException("Incorrect number of ensembles")
    if not 'mrbump_results' in ad or not len(ad['mrbump_results']): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    if not ad['mrbump_results'][0]['SHELXE_CC'] > 25: raise AmpleException("SHELXE_CC criteria not met")
    return
        
test_dict['rosetta_modelling'] = { 'args' : args_rosetta_modelling,
                                   'test' :  test_rosetta_modelling }

###############################################################################
#
# test from pre-existing models
#
###############################################################################
args_from_existing_models = args_vanilla + [
    '-models', '../../tests/testfiles/models',
]

def test_from_existing_models(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ad['ensembles'] or not len(ad['ensembles']) == 12: raise AmpleException("Incorrect number of ensembles")
    if not ('mrbump_results' in ad or len(ad['mrbump_results'])): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    if not ad['mrbump_results'][0]['SHELXE_CC'] > 25: raise AmpleException("SHELXE_CC criteria not met")
    return
        
test_dict['from_existing_models'] = { 'args' : args_from_existing_models,
                                     'test' :  test_from_existing_models }

###############################################################################
#
# test from quark models (also used as an opportunity to test the benchmark mode)
#
###############################################################################
args_from_quark_models = args_vanilla + [
    '-models', '../../tests/testfiles/decoys_200.tar.gz',
    '-native_pdb', '1DTX.pdb'                                         
]
 
def test_from_quark_models(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ad['ensembles'] or not len(ad['ensembles']) == 18: raise AmpleException("Incorrect number of ensembles")
    if not ('mrbump_results' in ad or len(ad['mrbump_results'])): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    return
         
test_dict['from_quark_models'] = { 'args' : args_from_quark_models,
                                   'test' : test_from_quark_models }

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
