#!/bin/bash

# 
# Script to run ideal helix test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

$CCP4/bin/ample \
-fasta input/2OVC.fasta \
-mtz input/2OVC-cad.mtz \
-use_shelxe True \
-ideal_helices True \
-nproc 8 \
-show_gui True \
 
