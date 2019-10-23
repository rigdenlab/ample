#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
#export PATH=\
#/opt/shelx:\
#$PATH

$CCP4/bin/ample \
-fasta input/3DCY.fasta \
-mtz input/3dcy-sf.mtz \
-homologs True \
-models input  \
-nproc 8 \
-show_gui True \

# Add below for running in benchmark mode
#-mustang_exe /opt/MUSTANG_v3.2.2/bin/mustang-3.2.1 \
#-alignment_file input/testthree.afasta \
#-native_pdb input/1DCY.pdb \
