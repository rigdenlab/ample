#!/bin/bash

# 
# Script to run toxd test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
ROSETTA_DIR=/opt/rosetta-3.5

export PATH=\
/opt/spicker:\
/opt/theseus_src:\
/opt/scwrl4:\
/opt/maxcluster:\
/opt/shelx:\
$PATH

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=/opt/ample-dev1

python ${AMPLEDIR}/bin/ample \
      -rosetta_dir $ROSETTA_DIR \
      -fasta ${PWD}/toxd_.fasta \
      -mtz ${PWD}/1dtx.mtz \
      -frags_3mers ${PWD}/aat000_03_05.200_v1_3 \
      -frags_9mers ${PWD}/aat000_09_05.200_v1_3 \
      -nmodels 30 \
      -percent 50 \
      -use_shelxe True \
      -shelx_cycles 10 \
      -use_arpwarp False \
      -use_buccaneer False \
      -molrep_only True \
      -nproc 2 \

      #-name TOXD \
      #-models_dir ${AMPLEDIR}/examples/toxd-example/ROSETTA_MR_0/models \
      #-ensembles_dir ${AMPLEDIR}/examples/toxd-example/ROSETTA_MR_0/ensembles_1 \
