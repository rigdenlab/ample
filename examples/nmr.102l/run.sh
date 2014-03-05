#!/bin/sh

# NMR ensembling example - 2LC9 is an ensemble model of a minor and transiently formed state of a T4 lysozyme mutant. The target 102l is X-ray data

ccp4-python -u ${CCP4}/bin/ample \
      -mtz 102l.mtz \
      -F F -SIGF SIGF -FREE FreeR_flag \
      -fasta 102L.fasta \
      -name 102l \
      -NMR_model_in 2LC9.pdb \
      -NMR_Truncate_only True \
      -use_arpwarp False \
      -molrep_only True \
      -nproc 1
