#!/usr/bin/env python

#

import signal, time
import subprocess
import os
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen
import glob, re


def get_shel_score(pdb):
  print 'here'
  pdb=open(pdb)
  score = 0

  pattern=re.compile('CC\s*=\s*(\d*\.\d*)')
  for line in pdb:
  
   if re.search(pattern, line):
     
     split = re.split(pattern,line)
     
     score = split[1]
  return score
# # # # # # # # # # # # # # # # 

def rank_results(mrbump_path, overpath):
 
 order = []

 if not os.path.exists( os.path.join(overpath, 'RESULTS')): 
   os.mkdir(os.path.join(overpath, 'RESULTS'))

 

 for folder in os.listdir(mrbump_path):
  if  os.path.isdir(mrbump_path+'/'+folder):
    

    if os.path.exists(mrbump_path+'/'+folder+'/phaser_shelx'): # phaser shelx
       if os.path.exists(mrbump_path+'/'+folder+'/phaser_shelx/orig.pdb'):
          score = get_shel_score(mrbump_path+'/'+folder+'/phaser_shelx/orig.pdb')       
          order.append([mrbump_path+'/'+folder+'/phaser_shelx/XYZOUT', float( score), mrbump_path+'/'+folder+'/phaser_shelx/HKLOUT' ])

    if os.path.exists(mrbump_path+'/'+folder+'/molrep_shelx'): # phaser shelx
       if os.path.exists(mrbump_path+'/'+folder+'/molrep_shelx/orig.pdb'):
          score = get_shel_score(mrbump_path+'/'+folder+'/molrep_shelx/orig.pdb')
          order.append([mrbump_path+'/'+folder+'/molrep_shelx/XYZOUT', float( score), mrbump_path+'/'+folder+'/molrep_shelx/HKOUT' ])
          print 'GOT'

 

 index=1
 order.sort(key=lambda tup: tup[1])
 order.reverse()
 for x in order:
   print x
   os.system('cp  ' + x[0] + ' '+ os.path.join(overpath, 'RESULTS', 'result_'+str(index)+'.pdb') )
   os.system('cp  ' + x[2] + ' '+ os.path.join(overpath, 'RESULTS', 'result_'+str(index)+'.mtz') )
 
   index+=1


if __name__ == "__main__":
 
 mrbump_path = '/home/rmk65/projects/hammer/test/ROSETTA_MR_3/MRBUMP'
 overpath = '/home/rmk65/projects/hammer/test/ROSETTA_MR_3'
 rank_results(mrbump_path, overpath)
