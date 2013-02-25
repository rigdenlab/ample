#!/bin/bash

# 
# Script to run toxd test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
ROSETTA_DIR=/opt/rosetta3.4

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
export PYTHONPATH=${AMPLEDIR}/python:$PYTHONPATH

${AMPLEDIR}/bin/ample \
      -fasta ${PWD}/toxd_.fasta \
      -MTZ ${PWD}/1dtx.mtz \
      -name TOXD \
      -frags3mers ${PWD}/aat000_03_05.200_v1_3 \
      -frags9mers ${PWD}/aat000_09_05.200_v1_3 \
      -rosetta_dir $ROSETTA_DIR \
      -run_dir $PWD \
      -make_frags False \
      -CLUSTER False \
      -NumShelxCyc 10 \
      -nmodels 30 \
      -TRYALL False  \
      -percent 50 \
      -use_shelxe True \
      -arpwarp False \
      -Buccaneer False \
      -molreponly True \
      -nproc 1 \
      -allatom True \
      -noClusters 1

#      -MakeModels False \
#      -F FP -SIGF SIGFP -FREE FreeR_flag \
