#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

#$CCP4/bin/ample.py \
/opt/ample-dev1/bin/ample.py \
-fasta 3DCY.fasta \
-mtz 3dcy-sf.mtz \
-homologs True \
-mustang_exe /opt/MUSTANG_v3.2.2/bin/mustang-3.2.1 \
-models . \
-nproc 8 \

# Add below for running in benchmark mode
#-alignment_file testthree.afasta \
#-native_pdb  1DCY.pdb \
