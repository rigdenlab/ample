#!/bin/bash

# 
# Script to run NMR remodelling test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

# Need to say where we can find rosetta
rosetta_dir=/opt/rosetta-3.5

$CCP4/bin/ample \
-name 1t00 \
-mtz input/1t00.mtz \
-fasta input/1T00.fasta  \
-rosetta_dir $rosetta_dir \
-nmr_model_in input/2DIZ.pdb \
-nmr_remodel True \
-frags_3mers input/1t00.200.3mers \
-frags_9mers input/1t00.200.9mers \
-nmr_process 1 \
-quick_mode True \
-nproc 4 \

