#!/usr/bin/python2.6

#

import signal, time
import subprocess
import os,re, sys


def get_rfree(pdb):
  pdb = open(pdb)

  freer_pattern = re.compile('REMARK\s*\d*\s*FREE R VALUE\s*\:\s*(\d*\.\d*)')
  FreeR = 'nan' 

  for line in pdb:

    if re.search(freer_pattern, line):
    #  print line
      split=  re.split(freer_pattern, line)
   #   print split
      FreeR = split[1]

  return FreeR
###########################
def get_shelx_score(RESULT):

    Best_result = 0
    result = open(RESULT)
    fail_pattern = ('\*\* Unable to trace map - giving up \*\*')
    pattern = re.compile('CC for partial structure against native data =\s*(\d*\.\d*)\s*%')
    for line in result: 
       if re.search(pattern, line):
        split= re.split(pattern, line)
      #  print split
        if float(split[1]) > Best_result:
          Best_result = float(split[1])
       if re.search(fail_pattern, line):
        Best_result = line.rstrip('\n')

    return Best_result
  
##########################
def make_log(mr_bump_path, final_log):
  final_log = open(final_log, "w")

  list_dir = os.listdir(mr_bump_path)
  for a_dir in list_dir:
    if os.path.isdir(mr_bump_path + '/'+a_dir):

     name=re.sub('search_', '', a_dir)
     name=re.sub('_mrbump', '', name)
     
     phaser_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/pdbclip/refine/phaser/refmac_phaser_loc0_ALL_'+name+'_PDBCLP.pdb'
     molrep_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/pdbclip/refine/molrep/refmac_molrep_loc0_ALL_'+name+'_PDBCLP.pdb'
     phaser_log = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/pdbclip/mr/phaser/phaser_loc0_ALL_'+name+'_PDBCLP.1.pdb'

     shelx_phaser = mr_bump_path + '/'+a_dir+'/phaser_shelx/RESULT'
     shelx_molrep = mr_bump_path + '/'+a_dir+'/molrep_shelx/RESULT'

     if os.path.exists(phaser_pdb):
      phaser_FreeR =   get_rfree(phaser_pdb)  
      phaser_shelx =   get_shelx_score(shelx_phaser)
      final_log.write('Ensembe ' + name+'  phaser FreeR: '+phaser_FreeR +  '  shelx score: '+str(phaser_shelx) + '\n')
      final_log.flush()
  
     if os.path.exists(molrep_pdb):
      molrep_FreeR =   get_rfree(molrep_pdb) 
      molrep_shelx =   get_shelx_score(shelx_molrep)  
      final_log.write('Ensembe ' + name+'  molrep FreeR: '+molrep_FreeR +  '  shelx score: '+str(molrep_shelx) + '\n')
      final_log.flush()
     

   




############
#mr_bump_path = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP'
#final_log = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP/FINAL'
#make_log(mr_bump_path,final_log)





