#!/bin/bash

# 
# Script to run missing domains test case
#

# Specify where rosetta is located
rosetta_dir=/opt/rosetta3.5

# add shelxe executable to path
export PATH=\
/opt/shelx:\
$PATH


$CCP4/bin/ccp4-python -m ample \
-mtz input/1k04_cad-unique.mtz \
-fasta input/1k04_.fasta  \
-domain_all_chains_pdb input/Known_40.pdb \
-missing_domain True \
-frags_3mers input/aa1k04_03_05.200_v1_3 \
-frags_9mers input/aa1k04_09_05.200_v1_3 \
-rosetta_dir $rosetta_dir \
-nproc 5 \


