#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
#export PATH=\
#/opt/shelx:\
#$PATH

# Path to the rosetta directory
rosetta_dir=/opt/rosetta_bin_linux_2015.39.58186_bundle/

$CCP4/bin/ample \
-fasta input/toxd_.fasta \
-models  ../../testfiles/decoys.tar.gz \
-mtz input/1dtx.mtz \
-nmodels 30 \
-percent 50 \
-num_clusters 1 \
-use_shelxe True \
-nproc 5 \
-show_gui True \

# Additional optional flags
# Add below to run from pre-made models
#-models ../../testfiles/models \
#-rosetta_dir $rosetta_dir \
#-frags_3mers input/aat000_03_05.200_v1_3 \
#-frags_9mers input/aat000_09_05.200_v1_3 \

# Thes are QUARK models

# Add below for running in benchmark mode
#-native_pdb  inpuy/1DTX.pdb \

# Add below for running from pre-made ensembles
#-ensembles ./ROSETTA_MR_0/ensembles_1 \

