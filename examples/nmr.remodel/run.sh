#!/bin/bash
# 
# Script to run NMR test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
rosetta_dir=/opt/rosetta-3.5
rosetta_dir=/data2/opt/rosetta_bin_linux_2015.05.57576_bundle

export PATH=\
/opt/maxcluster:\
/opt/shelx:\
$PATH

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=$CCP4/share/ample
AMPLEDIR=/home/jmht/ample-dev1

${AMPLEDIR}/bin/ample.py \
      -name 1t00 \
      -mtz 1t00.mtz \
      -fasta 1T00.fasta  \
      -rosetta_dir $rosetta_dir \
      -nmr_model_in 2DIZ.pdb \
      -nmr_remodel True \
      -frags_3mers 1t00.200.3mers \
      -frags_9mers 1t00.200.9mers \
      -nproc 8 \
      -submit_cluster True \
      -submit_qtype sge \

#      -quick_mode True \
#      -nmr_process 2 \
