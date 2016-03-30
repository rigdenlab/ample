#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

# Path to the rosetta directory
rosetta_dir=/opt/rosetta-3.5

$CCP4/bin/ccp4-python -m ample \
-rosetta_dir $rosetta_dir \
-fasta input/toxd_.fasta \
-mtz input/1dtx.mtz \
-frags_3mers input/aat000_03_05.200_v1_3 \
-frags_9mers input/aat000_09_05.200_v1_3 \
-nmodels 30 \
-percent 50 \
-use_shelxe True \
-nproc 5 \

# Additional optional flags
# Add below to run from pre-made models
#-models ../../tests/testfiles/models \

# Thes are QUARK models
#-models  ../../tests/testfiles/decoys_200.tar.gz \

# Add below for running in benchmark mode
#-native_pdb  inpuy/1DTX.pdb \

# Add below for running from pre-made ensembles
#-ensembles ./ROSETTA_MR_0/ensembles_1 \

