#!/bin/bash
# 
# Script to run NMR test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
rosetta_dir=/opt/rosetta-3.5

export PATH=\
/opt/shelx:\
$PATH

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=$CCP4/share/ample

${AMPLEDIR}/bin/ample.py \
      -name 1t00 \
      -mtz 1t00.mtz \
      -fasta 1T00.fasta  \
      -rosetta_dir $rosetta_dir \
      -nmr_model_in 2DIZ.pdb \
      -nmr_remodel True \
      -frags_3mers 1t00.200.3mers \
      -frags_9mers 1t00.200.9mers \
      -nmr_process 1 \
      -quick_mode True \
      -nproc 4 \

