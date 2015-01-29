#!/bin/bash

# 
# Script to run toxd test case
#

# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
export PATH=\
/opt/mrbump-trunk/bin:\
/opt/shelx:\
$PATH

# Path to the rosetta directory
rosetta_dir=/opt/rosetta-3.5

#$CCP4/bin/ample.py \
/opt/ample-dev1.testset/bin/ample.py \
-rosetta_dir $rosetta_dir \
-fasta toxd_.fasta \
-mtz 1dtx.mtz \
-frags_3mers aat000_03_05.200_v1_3 \
-frags_9mers aat000_09_05.200_v1_3 \
-num_clusters 1 \
-nmodels 30 \
-percent 50 \
-quick_mode True \
-use_shelxe True \
-shelxe_rebuild False \
-nproc 2 \
-models_dir ./models \

#-models_dir ./models \
#-ensembles_dir ./ROSETTA_MR_0/ensembles_1 \


