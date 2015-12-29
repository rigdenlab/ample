#!/usr/bin/env ccp4-python
'''
Created on 29 Dec 2015

@author: jmht
'''
import os, sys

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
    return script


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
    '-nproc', '2',
    '-no_gui', 'True',
    #'-do_mr', 'False',
]

args_rosetta = args_vanilla + [
    '-rosetta_dir', '/opt/rosetta-3.5',
    '-frags_3mers', 'aat000_03_05.200_v1_3',
    '-frags_9mers', 'aat000_09_05.200_v1_3',
    '-nmodels', '30',
]

# test from pre-existing models
args_models = args_vanilla + [
    '-models', '../../tests/testfiles/models',
]

# Test from quark models (also used as an opportunity to test the benchmark mode
args_quark = args_vanilla + [
    '-models', '../../tests/testfiles/decoys_200.tar.gz',
    '-native_pdb', '1DTX.pdb'
]

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