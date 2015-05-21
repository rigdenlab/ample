#!/bin/bash

# 
# Script to run missing domains test case
#

# Specify where rosetta is located
rosetta_dir=/opt/rosetta3.4

# add shelxe executable to path
export PATH=\
/opt/shelx:\
$PATH


${CCP4}/bin/ample.py \
-mtz 1k04_cad-unique.mtz \
-fasta 1k04_.fasta  \
-domain_all_chains_pdb  Known_40.pdb \
-missing_domain True \
-frags3mers aa1k04_03_05.200_v1_3 \
-frags9mers aa1k04_09_05.200_v1_3 \
-rosetta_dir $rosetta_dir \
-nproc 5 \


