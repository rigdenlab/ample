#!/bin/bash

$CCP4/bin/ample.py \
-fasta 1ujb.fasta \
-mtz 1ujb-sf.mtz \
-single_model 3c7t.pdb \
-do_mr False \
-truncation_method scores \
-truncation_scorefile 3c7t_scores.csv \
-truncation_scorefile_header residue Concoord AMN  \

