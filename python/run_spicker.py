#!/usr/bin/env python

#

import re
import os, glob
import sys
import subprocess
import time
import operator
import shutil


################
def check_pid_by_kill(pid):        
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True
##################
def check_pid(pid):    
 if os.path.exists('/proc/'+str(pid)):
   return True
 else:
   return False

###############
def get_length(pdb):
  pdb = open(pdb)
  counter = 0
  for line in pdb:
     pdb_pattern = re.compile('^ATOM')
     pdb_result = pdb_pattern.match(line)
     if pdb_result:
          atom = line[13:16]
          if re.search('CA', atom):
             counter+=1

 # print counter
  return str(counter)

#################
def run_spicker(path, outpath):

 """
 jmht
 Create the input files required to run spicker
 
 path: PATH_TO_MODELS
 outpath: RunDir+'/spicker_run'
 
 (See notes in spicker.f FORTRAN file for a description of the required files)
 
 """
 
  
 curr_dir = os.getcwd()
 working_dir = outpath 
 os.system('mkdir ' + working_dir)

 os.chdir(working_dir)
 pdbname = ''
 list_string = ''
 counter = 0

 # Input file for spicker with coordinates of the CA atoms for each of the PDB structures
 read_out = open(working_dir + '/rep1.tra1', "w")
 
 #jmht - a list of the full path of all PDBs - used so we can loop through it and copy the selected
 # ones to the relevant directory after we have run spicker
 file_list = open(working_dir + '/file_list', "w")
 
 for infile in glob.glob( os.path.join(path,  '*.pdb') ):
  name = re.split('/', infile)
  
  pdbname = str(name.pop())
  file_list.write(infile + '\n')
  list_string = list_string + pdbname+ '\n'
  counter +=1
  read = open(infile)
  
  length = get_length(infile)
  
  read_out.write('\t' + length + '\t926.917       '+str(counter)+'       '+str(counter)+'\n')
  for line in read:
     #print line
     pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)')
     result = re.match(pattern, line)
     if result:
       split = re.split(pattern, line)
       #print line
       #print split
      # print '\t' + split[5] + '\t' + split[6] + '\t' +split[7] 
       if split[2] == 'CA':
        read_out.write( '     ' + split[6] + '     ' + split[7] + '     ' +split[8] + '\n' )

  read.close()       
 file_list.close()
 #print list_string
#make rmsinp
 #length = get_length(path + '/' + pdbname)

 
 rmsinp = open('rmsinp', "w")
 rmsinp.write('1  ' + length + '\n\n')
 rmsinp.write(length + '\n')
 rmsinp.close()
#make tra.in
 tra = open('tra.in', "w")
 tra.write('1 -1 1 \nrep1.tra1')
           
 tra.close()    

 # Create the file with the sequence of the PDB structures
 seq = open('seq.dat', "w")
 a_pdb = open(path + '/' + pdbname)
 for line in a_pdb:
  #print line
  pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(\d*)\s')
  result = re.match(pattern, line)
  if result:
     split = re.split(pattern, line)
     #print line
    # print split
     if split[2] == 'CA':
      seq.write('\t' +split[5] + '\t' + split[3] + '\n')
################

def read_log_near_centroids(clus):
     """
     only get those closest to centroid
     jmht
     Select either the first 200 or all if there are less that 200
     """
     index = []
     cur_dir = os.getcwd()
     log = open(cur_dir+'/str.txt')
     pattern = re.compile('\s*(\d*)\s*.*\s*(\d*)\s*rep1.tra1')
     for line in log:
       if re.match(pattern, line):
         split = re.split('\s*', line)
         if split[1] != '':
          if len(split)>6:
    
           if int(split[1]) == clus:
            #put into touples
    
    
            index.append( [int(split[6]) , float(split[4])] )
           # print line
           # print split
            
    
     #filter and retrun
     KEEP  = min(200, len(index))  #keep the first 200, or all
    
     closest=[]
     K=0
     sorted_by_second = sorted(index, key=lambda tup: tup[1])
     while K<KEEP:
      # print sorted_by_second[K]
       closest.append(sorted_by_second[K][0])
       K+=1
    
    # print closest
    
     return closest

################
def read_log():
  index = []
  cur_dir = os.getcwd()
  log = open(cur_dir+'/rst.dat')
  pattern = re.compile('\s*(\d*)\s*.*\s*(\d*)\s*rep1.tra1')
  for line in log:
    if re.search('rep1.tra1', line):
         line = line.strip()
         s=re.split('(\s*)',line )
         if len(s)>2: 
          # print line
          # print s

           print 'cluster '+s[0]  +' has '+s[2] +' models'

  return 
##############################


def RUN_SPICKER(models, spicker_runpath, spickerexe, no_clusters_sampled, overpath):
   run_spicker(models, spicker_runpath)
   p = subprocess.Popen(spickerexe, shell=True, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.PIPE)
   p.wait()


   os.chdir(spicker_runpath)



   clus=1
   while clus < no_clusters_sampled+1:

     list_of_models = read_log_near_centroids(clus)  ### try only near centroid # use read_log for all
     models_paths = []
     read_log()
     if not os.path.exists(os.path.join(overpath, "S_clusters" )):
      os.mkdir(os.path.join(overpath, "S_clusters" ))
    
     if not os.path.exists( os.path.join(overpath, "S_clusters","cluster_"+str(clus) )):
      os.mkdir(os.path.join(overpath, "S_clusters","cluster_"+str(clus) ))

   #  print list_of_models
     i = 1
     models_list = open(spicker_runpath+'/file_list')
     for line in models_list:
       
        if i in list_of_models:
        # print line
         models_paths.append(line.rstrip('\n'))
        i+=1
     models_list.close()
    # print models_paths
     for eachpdb in models_paths:      
       shutil.copy(eachpdb, os.path.join(overpath, "S_clusters","cluster_"+str(clus) ))
     clus+=1




if __name__ == "__main__":
  
  models='/data2/jac45/nmr/ROSETTA_MR_3/models'
 # models='/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/VERIFY/Verify_spicker/SPICKER_script/models'
  spicker_run_dir='/data2/jac45/nmr/ROSETTA_MR_3/spicker_run'
  spickerexe='spicker'
  no_clusters_sampled=1
  #rosetta RunDir:
  overpath='/data2/jac45/nmr/ROSETTA_MR_3'
  

  RUN_SPICKER(models, spicker_run_dir, spickerexe, no_clusters_sampled, overpath)








