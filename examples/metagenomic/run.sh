#!/bin/bash

# 
# Script to run coiled-coil test case
#

# Set path to include where shelxe is located
#export PATH=\
#/opt/shelx:\
#$PATH

# Path to the rosetta directory
rosetta_dir=/opt/rosetta_bin_linux_2015.39.58186_bundle/

$CCP4/bin/ample \
-fasta input/5edl.fasta \
-mtz input/5edl.mtz \
-models ../../../testfiles/metagenomic_models \
-nproc 5 \
-show_gui True \

# Add below for running in benchmark mode
#-native_pdb  input/5edl.pdb \
