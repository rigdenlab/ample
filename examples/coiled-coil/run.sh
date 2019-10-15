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
-rosetta_dir $rosetta_dir \
-fasta input/1g1j.fasta \
-mtz input/1g1j.mtz \
-frags_3mers input/aa4u5t_03_05.200_v1_3 \
-frags_9mers input/aat000_09_05.200_v1_3 \
-multimer_modelling tetramer \
-nmodels 30 \
-nmasu 2 \
-shelxe_rebuild_buccaneer True \
-max_shexle_resolution 4 \
-nproc 5 \

# Additional optional flags
# Add below to run from pre-made models
#-models ../../testfiles/models \

# Add below for running in benchmark mode
#-native_pdb  input/1g1j.pdb \

# Add below for running from pre-made ensembles
#-ensembles ./ROSETTA_MR_0/ensembles_1 \

