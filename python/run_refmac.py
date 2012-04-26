#!/usr/bin/python2.6

#

import re
import glob
import sys, string, os 
########
def parse_file(filename):
    """
    Parse a file, returning a list of tags.
    Returns None on error.
    """
    
    f = open(filename,'r')

    tags = []
    
    for line in f:
        tags.append( line )
        
    f.close()
    return tags
##########
def refmac(xyzin, hklin):

  cur_dir = os.getcwd()
  SCORE = 'nan'
  tags = parse_file(xyzin)
  if len(tags) == 0:

    return os.path.join(cur_dir, 'refined.mtz'), SCORE
  mr_bump = open(cur_dir + '/refmac_run', "w") 
  mr_bump.write('#!/bin/sh\n')
  mr_bump.write('refmac5 xyzin ' + xyzin + ' hklin ' + hklin + ' hklout refined.mtz xyzout refined.pdb <<eof\n')
  mr_bump.write('ncyc 1\n')
  mr_bump.write('END\n')
  mr_bump.write('eof\n')
  mr_bump.close()
 
  os.system('chmod uoga=wrx '+cur_dir + '/refmac_run')
  os.system(cur_dir + '/refmac_run >refmac_out')  
#  print 'Done'

  out = open(cur_dir + '/refmac_out')
  pattern = re.compile('R free\s*(\d*\.\d*)\s*(\d*\.\d*)')
  for line in out:
   if re.search(pattern, line):
     #print line
     split = re.split(pattern, line)
     #print split
     SCORE = split[2]
  
  return os.path.join(cur_dir, 'refined.mtz'), os.path.join(cur_dir, 'refined.pdb'), SCORE


#########
#xyzin = '/home/jaclyn/Rosetta_bucc/rosetta_warp/1EJG_1/molrep_ALL_ATOM_trunc1_Rad3residue_buccaneer1.pdb' 
#hklin = '/home/jaclyn/Rosetta_bucc/rosetta_warp/1EJG_1/refmac_phaser_HKLOUT_loc0_ALL_POLYALA_trunc1rad_ALIGNED_rad_1_PDBCLP.mtz' 

#refmac(xyzin, hklin)
