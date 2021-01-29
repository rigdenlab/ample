#!/bin/bash

$CCP4/bin/ample \
-fasta $PWD/input/1ujb.fasta \
-mtz $PWD/input/1ujb-sf.mtz \
-single_model $PWD/input/3c7t.pdb \
-do_mr False \
-truncation_method scores \
-truncation_scorefile $PWD/input/3c7t_scores.csv \
-truncation_scorefile_header residue Concoord AMN  \

