#!/usr/bin/python2.6

#

import re
import os, glob
import sys
import subprocess
import time
import operator
import argparse
import random
from multiprocessing import Process, Queue, JoinableQueue, Pool, Value, Array
import pickle
import copy

#
# This Adds All sicechains to the Pdb files
#

def add_sidechains_SCWRL(SCWRL,path, outpath, prefix ):
  Path_to_SCWRL = SCWRL

  string = ''
  infiles = []

  for infile in glob.glob( os.path.join(path, '*.pdb') ):
    string  = string + infile + ' '
    infiles.append(infile) 
    #print infile
########################################## ADD ALL

  for each_file in infiles:
   residue_list = []
   previous = 0
   name = re.split('/', each_file)
   pdbname = str(name.pop()) 
   my_infile = open (each_file)
   for pdbline in my_infile:
     pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
     pdb_result = pdb_pattern.match(pdbline)
     if pdb_result:
      pdb_result2 = re.split(pdb_pattern, pdbline )
      #print pdb_result2[5] + ' |' + pdb_result2[3] + '|'
      current = previous +1
      if int(pdb_result2[5]) == current:               ############# make sequence lise
       #print  str(previous) + ' next residue ' + pdb_result2[5]
       previous = previous +1
       if pdb_result2[3] == 'ARG':
         residue_list.append('R')     ### caps = add sidechains
       if pdb_result2[3] == 'GLN':
         residue_list.append('Q')
       if pdb_result2[3] == 'PHE':
         residue_list.append('F')
       if pdb_result2[3] == 'TYR':
         residue_list.append('Y')
       if pdb_result2[3] == 'TRP':
         residue_list.append('W')
       if pdb_result2[3] == 'HIS':
         residue_list.append('H')
       if pdb_result2[3] == 'LYS':
         residue_list.append('K')
       if pdb_result2[3] == 'GLY':
         residue_list.append('G')
       if pdb_result2[3] == 'ALA':
         residue_list.append('A')
       if pdb_result2[3] == 'SER':
         residue_list.append('S')
       if pdb_result2[3] == 'PRO':
         residue_list.append('P')
       if pdb_result2[3] == 'GLU':
         residue_list.append('E')
       if pdb_result2[3] == 'ASP':
         residue_list.append('D')
       if pdb_result2[3] == 'THR':
         residue_list.append('T')
       if pdb_result2[3] == 'CYS':
         residue_list.append('C')
       if pdb_result2[3] == 'MET':
         residue_list.append('M')
       if pdb_result2[3] == 'LEU':
         residue_list.append('L')
       if pdb_result2[3] == 'ASN':
         residue_list.append('N')
       if pdb_result2[3] == 'ILE':
         residue_list.append('I')
       if pdb_result2[3] == 'VAL':
         residue_list.append('V')
       #else:
        #print ' UNRECOGNISED AA ' + pdb_result2[3]
        #sys.exit()
#   print residue_list    
   res_string = ''
   for i in residue_list:
    res_string = res_string + i
   my_seqfile = open (path + 'seq', "w")
#   print 'string ' + res_string 
   my_seqfile.write(res_string)
   my_seqfile.close()

   cmd = Path_to_SCWRL + ' -i ' + each_file + ' -o ' + outpath + '/' + prefix+ '_'+pdbname + ' -s ' + path + 'seq' ### add sidechains
#   print cmd
   os.system (cmd )

   edit_file=open(outpath + '/' + prefix+ '_' +pdbname)    ### remove headers
   lines_list= edit_file.readlines()
   edit_file.close()

   change_file=open(outpath + '/' +prefix+ '_' + pdbname, "w")

   for line in lines_list:
#     print line
     pdb_pattern = re.compile('^ATOM')
     pdb_result = pdb_pattern.match(line)
     if pdb_result:
#      print 'match'
      change_file.write(line)

