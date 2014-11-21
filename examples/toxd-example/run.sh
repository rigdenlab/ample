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

#$CCP4/bin/ample.py \
/opt/ample-dev1/bin/ample.py \
-rosetta_dir /opt/rosetta-3.5 \
-fasta toxd_.fasta \
-mtz 1dtx.mtz \
-frags_3mers aat000_03_05.200_v1_3 \
-frags_9mers aat000_09_05.200_v1_3 \
-nmodels 30 \
-percent 50 \
-quick_mode True \
-use_shelxe True \
-nproc 6 \
-models_dir ./models \

#-models_dir ./models \
#-ensembles_dir ./ROSETTA_MR_0/ensembles_1 \


