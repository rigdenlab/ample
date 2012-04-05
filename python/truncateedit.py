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
import subprocess
import shlex
import shutil

import cluster_with_MAX
import SCWRL_edit
#############cluster_with_MAX####
def make_list_to_keep(theseus_out, THESEUS_threthold): #make a list of residues to keep under variance threshold
 add_list =[]

 theseus_out = open(theseus_out)
 for line in theseus_out:
  pattern = re.compile('^(\d*)\s*(\w*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
  result = pattern.match(line)
  if result:
   seq = re.split(pattern, line)
   #print seq
   if (seq[6]) != 'T':

    if float(seq[4]) <= float(THESEUS_threthold):
     if seq[1] != '':
    #  print 'GOT ' + str(seq[4]) + '  thresh ' + str(THESEUS_threthold) + ' KEEP ' + str(seq[3]) 
      add_list.append(seq[3])

 print add_list
 return add_list
################################
def Run_Rosetta(string, no_of_files, radius, Rosetta_cluster, RDB):  # rosetta can fail so repeat job
  rosettaout = os.getcwd()
  condition = 0
  cluster_logfile="collect_out.log"

  while condition == 0:
    #os.system (Rosetta_cluster+ ' -database '+RDB+' -cluster:radius ' +str(radius) + ' -in:file:fullatom '  +  ' -in:file:s '+string+' > collect_out')

    command_line = Rosetta_cluster+ ' -database '+RDB+' -cluster:radius ' +str(radius) + ' -in:file:fullatom ' + ' -in:file:s '+string 

    process_args = shlex.split(command_line)
    p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE)

    (child_stdout, child_stdin) = (p.stdout, p.stdin)

    # Write the keyword input
    child_stdin.close()

    # Open the log file for writing
    cluster_log=open(cluster_logfile, "a")

    # Watch the output for successful termination
    out=child_stdout.readline()

    if os.environ["ROSETTA_VERSION"] == "3.2.1":
       file_pattern = re.compile('Structures\:\s*(\d*)')
    if os.environ["ROSETTA_VERSION"] == "3.1":
       file_pattern = re.compile('^Clustering\s*(\d*)\s*structures')
    while out:
       if file_pattern.search(out):
          file_result2 = re.split(file_pattern, out)
          if int(file_result2[1]) !=  no_of_files:
             print 'fail there are ' + str(no_of_files) + ' files, got ' + file_result2[1]
          if int(file_result2[1]) ==  no_of_files:
             print ' SUCCESS there are ' + str(no_of_files) + ' files, got ' + file_result2[1] 
             condition = 1

       cluster_log.write(out)
       cluster_log.flush()
       out=child_stdout.readline()

    child_stdout.close()
    cluster_log.close()


    #my_sdtout = open (rosettaout + '/collect_out') 
    #for line in my_sdtout:
     #file_pattern = re.compile('^Clustering\s*(\d*)\s*structures')
    # file_pattern = re.compile('Structures\:\s*(\d*)')
    # file_result = file_pattern.match(line)
    # if file_result:
    #  print line
#
#      file_result2 = re.split(file_pattern, line )
#      print file_result2
#      if int(file_result2[1]) !=  no_of_files:
#       print 'fail there are ' + str(no_of_files) + ' files, got ' + file_result2[1]
#      if int(file_result2[1]) ==  no_of_files:
#       print ' SUCCESS there are ' + str(no_of_files) + ' files, got ' + file_result2[1] 
#       condition = 1
####################################
def try_theseus(cmd):#  theseus can fail so try Run a command with a timeout after which it will be forcibly killed.

 has_worked = False
 while has_worked == False:

   p = subprocess.Popen(cmd, shell = True)
   time.sleep(5)
   print p.poll()
   if p.poll() is None: #still running      
      p.kill()
      print 'timed out'
   else:
     print p.communicate()
     has_worked = True


  



#################################
def Align_rosetta_fine_clusters_with_theseus(Rad_path, THESEUS, ):  #if a cluster is present: align cluster  if no cluster is present: try to align most similar
  no_of_files=0

  string = '' 
  for infile in glob.glob( os.path.join(Rad_path, '*.pdb') ):  #check how in top cluster
    name = re.split('/', infile)
    pdbname = str(name.pop())
    #print pdbname
    cluster_name = re.split('\.', pdbname)

    if int(cluster_name[1]) == 0:
     if no_of_files < 30: ######################### LIMIT number in ensemble
      string  = string + infile + ' '
      no_of_files = no_of_files + 1

  print no_of_files
  if no_of_files >5:
     try_theseus(THESEUS + ' -r ' + 'cluster_' + str(0) + ' -a0 ' +string)

  #   if clustering has failed, try anyway
  no_of_files=0
  string = ''  
  if no_of_files <5:
    for infile in glob.glob( os.path.join(Rad_path, '*.pdb') ):
      no_of_files = no_of_files + 1
      if no_of_files < 30: ######################### LIMIT number in ensemble
        string  = string + infile + ' '
    try_theseus(THESEUS + ' -r ' + 'cluster_' + str(0) + ' -a0 ' +string)


################################
def make_ensembles(trunc_path, Rosetta_cluster, threshold,THESEUS, RDB):

  ensembles_made = []

  RADS = [1,2,3] # radius thresholds
  no_files = 0
  string = ''
  for infile in glob.glob( os.path.join(trunc_path, '*.pdb') ):
     string  = string + infile + ' '
     no_files +=1

  for RAD in RADS:
    Rad_path = trunc_path + '/fine_clusters_'+str(RAD)
    os.system('mkdir ' +trunc_path + '/fine_clusters_'+str(RAD))
    os.chdir(trunc_path + '/fine_clusters_'+str(RAD))
    Run_Rosetta(string, no_files,  RAD, Rosetta_cluster, RDB)
    
    ensemble_path = trunc_path + '/fine_clusters_'+str(RAD)+'_ensemble'
    os.system('mkdir '+ensemble_path )
    os.chdir(ensemble_path)
    Align_rosetta_fine_clusters_with_theseus(Rad_path, THESEUS, )
    os.system('mv ' +ensemble_path + '/cluster_0_sup.pdb '+ensemble_path + '/trunc_'+str(threshold)+'_rad_' +str(RAD)+'.pdb')
    ensembles_made.append(ensemble_path + '/trunc_'+str(threshold)+'_rad_' +str(RAD)+'.pdb')
  return ensembles_made




################################

def truncate(THESEUS, models_path, out_path, Rosetta_cluster, RDB ): #truncate models in one folder
  all_ensembles = []
  thresholds = [1,2,3,4,5,6] ## Trucnation thresholds
  run_dir = os.getcwd()
  if not os.path.exists(out_path):
   os.mkdir(out_path)
  os.chdir(out_path)
 
  number_of_models = 0
  list_of_pdbs = []
  string = ''

 # for infile in glob.glob( os.path.join(models_path, '*.pdb') ):
 #    number_of_models +=1
 #    string  = string + infile + ' '
 #    list_of_pdbs.append(infile)


  #--------------------------------
  # get variations between pdbs
  #--------------------------------  
  os.chdir(models_path)
  print run_dir

  os.system (THESEUS + ' -a0 `ls *.pdb | xargs`  >theseus_data')
  print 'done'
  #shutil.move('theseus_variances.txt', run_dir)
  os.system('mv theseus* '+run_dir)
  #time.sleep(1)
  print run_dir
  T_data = out_path + '/theseus_data'

  #-------------------------------
  #truncate
  #----------------------------------
  truncate_log = open(out_path + '/truncate_log', "w")
  truncate_log.write('This is the number of residues kept under each truncation threshold\n\nthreshold\tnumber of residues\n')
  for threshold in thresholds:
    T_data = out_path + '/theseus_variances.txt'
    add_list = make_list_to_keep(T_data, threshold)
    truncate_log.write( str(threshold) +'\t' + str(len(add_list)) + '\n' )
    trunc_out=''

    continue_trunc = True
    if len(add_list) < 1:  #######   limit to size of pdb to keep, if files are too variable, all residues are removed
       
       continue_trunc = False

 
    if continue_trunc == True:
     print '-----truncatung at '+str(threshold)+'------'
     trunc_out = out_path + '/trunc_files_' + str(threshold)
     #print trunc_out
     os.system('mkdir ' +trunc_out)
     for infile in glob.glob( os.path.join(models_path,  '*.pdb') ):

      name = re.split('/', infile)
      pdbname = str(name.pop())

      my_infile = open(infile, "r")
      pdbout = trunc_out+'/' + pdbname
      pdb_out = open (pdbout , "w")
      for pdbline in my_infile:
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)
        if pdb_result:
          pdb_result2 = re.split(pdb_pattern, pdbline )

          for i in add_list : #convert to ints to compare

            if int(pdb_result2[5]) == int(i):

               pdb_out.write(pdbline)
      
      pdb_out.close()
     print 'making ensembles', threshold
     made_ens = make_ensembles(trunc_out,  Rosetta_cluster, threshold, THESEUS, RDB)  ### MAKE ensemble for this trucnation level
     for a_ens in made_ens:
       all_ensembles.append(a_ens)
  #print 'got ensembles', all_ensembles 
  return all_ensembles
###################
#/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus 
#/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0 



 #            /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0
#THESEUS =  '/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus'
#models_path='/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/VERIFY/Funnels/test_truncation'#
#out_path = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/VERIFY/Funnels/test_truncation/out'
#Rosetta_cluster = '/home/jaclyn/programs/rosetta3.1_Bundles/rosetta_source/bin/cluster.linuxgccrelease' 
#RDB='/home/jaclyn/programs/rosetta3.1_Bundles/rosetta_database'
#list_of_ensembles = truncate(THESEUS, models_path, out_path,Rosetta_cluster, RDB  )



#for each_ens in list_of_ensembles:
# SCWRL_edit.edit_sidechains(each_ens, '/home/jaclyn/Baker/test_cases/1AAR/ensembles/')
