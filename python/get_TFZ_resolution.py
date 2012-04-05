#!/usr/bin/python2.6

#

import re
import os, glob
import sys
import subprocess
import time


def get_resolution(pdb):
 
  RESOLUTION = ''
  TFZ = ''
  FREE_R = ''

  pdb = open(pdb)
  for line in pdb: 

   get_res = re.compile('REMARK   3   RESOLUTION RANGE HIGH \(ANGSTROMS\) \:\s*(\d*.\d*)' )
   res_result  = get_res.search(line)
   if res_result:
    print line
    split_line = re.split(get_res, line)
    print split_line
    RESOLUTION = split_line[1]


   get_freer = re.compile('REMARK   3   FREE R VALUE                     :\s*(\d*.\d*)' )
   freer_result  = get_freer.search(line)
   if freer_result:
    print line
    split_line = re.split(get_freer, line)
    print split_line
    FREE_R = split_line[1]

  return RESOLUTION, FREE_R
##########################
def get_TFZ(pdb):

  TFZ = ''
  
  pdb = open(pdb)
  for line in pdb: 

   get_tfz = re.compile('REMARK  RFZ=(\d*.\d*)\s*TFZ=(\d*.\d*)' )
   tfz_result  = get_tfz.search(line)
   if tfz_result:
    print line
    split_line = re.split(get_tfz, line)
    print split_line
    TFZ = split_line[2]

  return TFZ

#tfz  = get_TFZ('/home/jaclyn/TASSER/RESULT/ITASSER_run/1OKS/search_1OKS_mrbump/data/loc0_ALL_1OKSA_3/pdbclip/mr/phaser/phaser_loc0_ALL_1OKSA_3_PDBCLP.1.pdb')

#print ' GOT ' + resolution + ' ' +free_r
