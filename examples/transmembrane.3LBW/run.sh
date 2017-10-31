#!/bin/bash

# 
# Script to run Transmembrane test case
#

# Need to specify where rosetta is located
rosetta_dir=/opt/rosetta_2014.35.57232_bundle

$CCP4/bin/ample \
-name 3LBW \
-mtz input/3LBW.sf.mtz   \
-fasta input/3LBW.fasta  \
-rosetta_dir $rosetta_dir \
-frags_3mers input/3LBW.200.3mers  \
-frags_9mers input/3LBW.200.9mers  \
-transmembrane True \
-contact_file input/c5313da6-ce36-47ea-81ee-79925b2fa836.metapsicov.stage1.txt \
-nmodels 20 \
-contact_format metapsicov \
-nproc 5 \
-do_mr False

