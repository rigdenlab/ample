#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''
import cPickle
import glob
import os, sys
import unittest

"""
Need:
name
arguments
test_function

To run:
nproc or submit_queue

each test needs to return the path to the resultsd.pkl

loop through directories
build up list of jobs and tests
write scripts, run and wait to complete
run tests
report results

scripts_and_tests = []
for d in dirs:
   scripts_and_tests.append(d)


for  
    

"""

sys.path.append('/opt/ample-dev1/python')
import workers

class AmpleException(Exception): pass

def write_script(name, args):
    ample = os.path.join(os.environ['CCP4'],'bin', 'ample.py')
    ample = os.path.join('/opt', 'ample-dev1', 'bin', 'ample.py')
    
    windoze = True if sys.platform.startswith('win') else False
    header = '' if windoze else '#!/bin/bash'
    ext = '.bat' if windoze else '.sh'
    
    script = name + ext
    with open(script, 'w') as f:
        f.write(header + os.linesep)
        f.write(os.linesep)
        f.write(ample + " \\" + os.linesep)
        # Assumption is all arguments are in pairs
        arg_list = [ " ".join(args[i:i+2]) for i in range(0, len(args), 2) ]
        f.write(" \\\n".join(arg_list))
        f.write(os.linesep)
        f.write(os.linesep)
    
    os.chmod(script, 0o777)
    return os.path.abspath(script)


# Go through each test directory
# 'homologs'
# 'ideal-helices'
# 'missing-domain.1k04'
# 'nmr.remodel'
# 'nmr.truncate'
# 'toxd-example'
# 'transmembrane.3LBW'


# vanilla test
args_vanilla = [
    '-fasta', 'toxd_.fasta',
    '-mtz', '1dtx.mtz',
    '-percent', '50',
    '-nproc', '1',
    '-no_gui', 'True',
    #'-do_mr', 'False',
]



testd = {}

###############################################################################
#
# Making models
#
###############################################################################
#/opt/rosetta_2014.35.57232_bundle
args_rosetta_modelling = args_vanilla + [
    '-rosetta_dir', '/opt/rosetta-3.5',
    '-frags_3mers', 'aat000_03_05.200_v1_3',
    '-frags_9mers', 'aat000_09_05.200_v1_3',
    '-nmodels', '30',
]


def test_rosetta_modelling(resultsd_pkl):
    with open(resultsd_pkl) as f: ad = cPickle.load(f)
    if not ad['models_dir'] or not len(glob.glob(ad['models_dir']+"/*.pdb'")) == 30: raise AmpleException("Incorrect number of models")
    if not ('ensembles' in ad or len(ad['ensembles'])): raise AmpleException("Incorrect number of ensembles")
    if not ('mrbump_results' in ad or len(ad['mrbump_results'])): raise AmpleException("No MRBUMP results")
    if not ad['success']: raise AmpleException("Job did no succeed")
    return
        
testd['rosetta_modelling'] = { 'args' : args_rosetta_modelling,
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
    return
        
testd['from_existing_models'] = { 'args' : args_from_existing_models,
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
        
testd['from_quark_models'] = { 'args' : args_from_quark_models,
                               'test' : test_from_quark_models }

###############################################################################
#
# End Test Setup
#
###############################################################################


# Create scripts and path to resultsd
scripts = []
for name in testd.keys():
    script = write_script(name, testd[name]['args'] + ['-work_dir', name])
    scripts.append(script)
    testd[name]['resultsd'] = os.path.join(name,'resultsd.pkl')

# Run all the jobs
nproc = 5
submit_cluster = False
submit_qtype = 'sge'
submit_array = True
if True:
    workers.run_scripts(job_scripts=scripts,
                        monitor=None,
                        chdir=True,
                        nproc=nproc,
                        job_time=3600,
                        job_name='test',
                        submit_cluster=submit_cluster,
                        submit_qtype=submit_qtype,
                        submit_queue=None,
                        submit_array=submit_array,
                        submit_max_array=None)


# Now run the tests
for name in testd.keys():
    try:
        testd[name]['test'](testd[name]['resultsd'])
        print "Job \'{0}\' succeeded".format(name)
    except AmpleException as ae:
        print "* Job \'{0}\' failed a test: {1}".format(name, ae)
    except Exception as e:
        print "*** Job \'{0}\' generated an exeption: {1}".format(name, e)

sys.exit()

write_script('making_rosetta_models', args_rosetta)
write_script('from_existing_models', args_models)
write_script('from_quark_models', args_models)


sys.path.insert(0,'/opt/ample-dev1/bin')

from ample import Ample

AMPLE = Ample()
#main(args_vanilla)
AMPLE.run(args_quark)

print "GOT SUCCESS ", AMPLE.amopt.d['success']

#-restart_pkl AMPLE_4/resultsd.pkl \
#-mrbump_dir AMPLE_0/MRBUMP \


# Add below for running from pre-made ensembles
#-ensembles ./ROSETTA_MR_0/ensembles_1 \

# Add below for running with contact predictions
#-contact_file toxd_.pconsc2.CASPRR \
#-contact_file toxd_.gremlin.CASPRR \
#-bbcontacts_file toxd_.bbcontacts.CASPRR \
#-constraints_file toxd_.cst \
#-native_pdb 1DTX_std.pdb \
#-energy_function FADE_default \
