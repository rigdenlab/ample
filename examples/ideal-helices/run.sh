#!/bin/bash

# 
# Script to run toxd test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

/opt/ample-dev1/bin/ample.py \
-fasta 2OVC.fasta \
-mtz 2OVC-cad.mtz \
-use_shelxe True \
-ideal_helices True \
-nproc 8 \
