#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

$CCP4/bin/ccp4-python -m ample \
-fasta input/2OVC.fasta \
-mtz input/2OVC-cad.mtz \
-use_shelxe True \
-ideal_helices True \
-nproc 8 \
