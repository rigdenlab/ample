#!/bin/sh

# NMR ensembling example - 2LC9 is an ensemble model of a minor and transiently formed state of 
# a T4 lysozyme mutant. The target 102l is X-ray data

# Set path to include where shelxe is located
export PATH=\
/opt/shelxe:\
$PATH

$CCP4/bin/ccp4-python -m ample \
-mtz input/102l.mtz \
-fasta input/102L.fasta \
-name 102l \
-nmr_model_in input/2LC9.pdb \
-quick_mode True

