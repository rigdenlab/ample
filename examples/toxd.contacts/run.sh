#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
#export PATH=\
#/opt/shelx:\
#$PATH

# Path to the rosetta directory
#rosetta_dir=/opt/rosetta-3.5
rosetta_dir=/opt/rosetta_src_2015.39.58186_bundle

#$CCP4/bin/ample.py \
../../bin/ample.py \
-rosetta_dir $rosetta_dir \
-fasta toxd_.fasta \
-mtz 1dtx.mtz \
-frags_3mers aat000_03_05.200_v1_3 \
-frags_9mers aat000_09_05.200_v1_3 \
-nmodels 5 \
-percent 50 \
-use_shelxe True \
-nproc 5 \
-contact_file toxd_.pconsc2.CASPRR \
-native_pdb  1DTX_std.pdb \
-no_gui True \
-do_mr False \

# Add below for running in benchmark mode
#-native_pdb  1DTX_std.pdb \

# Add below for running with contact predictions
#-contact_file toxd_.gremlin.CASPRR \
#-bbcontacts_file toxd_.bbcontacts.CASPRR \
#-constraints_file toxd_.cst \
#-native_pdb 1DTX_std.pdb \
#-energy_function FADE_default \
