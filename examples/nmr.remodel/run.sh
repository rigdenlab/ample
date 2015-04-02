#!/bin/bash
# 
# Script to run NMR test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
ROSETTA_DIR=/opt/rosetta_2014.35.57232_bundle

export PATH=\
/opt/maxcluster:\
/opt/shelx:\
$PATH

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=$CCP4/share/ample
AMPLEDIR=/opt/ample-dev1

${AMPLEDIR}/bin/ample.py \
      -name 1t00 \
      -mtz 1t00.mtz \
      -fasta 1T00_.fasta  \
      -rosetta_dir $ROSETTA_DIR \
      -NMR_model_in 1OKD.pdb \
      -NMR_remodel True \
      -frags3mers 1t00_.200.3mers \
      -frags9mers 1t00_.200.9mers \
      -quick_mode True \
      -nproc 5 \

