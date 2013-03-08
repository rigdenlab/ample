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
from subprocess import PIPE

import cluster_with_MAX

def One_model_only(list_of_ensembles, rundir):
 if not os.path.exists(rundir + '/Top_model_ensembles'):
   os.mkdir(rundir + '/Top_model_ensembles')  

 
 for pdb in list_of_ensembles:
  if os.path.exists(pdb):
   name = os.path.split(pdb)
   #print name
   #print pdb
   outpdb = open(rundir + '/Top_model_ensembles/'+name[-1], "w")

   for line in open(pdb):
       if re.search('ENDMDL', line):
           done = True
           break
       outpdb.write(line)
   outpdb.close()

 outlist = []
 for outpdb in os.listdir(rundir + '/Top_model_ensembles'):
     outlist.append( rundir + '/Top_model_ensembles/'+ outpdb)
 return  outlist

#############cluster_with_MAX####
def make_list_to_keep(theseus_out, THESEUS_threthold):
    """
    Make a list of residues to keep under variance threshold
    INPUTS:
    theseus_out: theseus_variances.txt output file
    THESEUS_threthold: the threshold variance
    
    RETURNS:
    a list of the residue indexes
    """
    add_list =[]
    
    theseus_out = open(theseus_out)
    for line in theseus_out:
     line = re.sub('RES\s*', '', line)
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
    
    #print add_list
    return add_list
#END make_list_to_keep

def chunks(a_list, percent_interval):
 

 for i in xrange(0, len(a_list), percent_interval):
        yield a_list[i:i+percent_interval ]
##############################
def fly_threshold(theseus_out, percent):
    """
    Make a list of residues to keep under variance threshold
    Threshold is chosen based on number of residues kept
    
    INPUT:
    theseus_out: theseus_variances.txt output file
    percent: 
    
    """
    add_list =[]
    
    # List of variances ordered by residue index
    var_list=[]
    
    # print theseus_out
    theseus_out = open(theseus_out)
    #jmht - build up a list of the variances of each residue
    for line in theseus_out:
     #print line
     line = re.sub('RES\s*', '', line)  #jmht - this probably not needed as we removed them earlier
     pattern = re.compile('^(\d*)\s*(\w*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
     result = pattern.match(line)
     if result:
     # print line
      #seq = re.split('\s*', line)
      seq = re.split(pattern, line)
     # print seq
     #jmht - I think this is just to skip the header line so should just do that at the start
      if not re.search('ATOM', line):
      #if (seq[6]) != 'T':
         #print seq[4]
    
         var_list.append(float(seq[4]))
    
    length = len(var_list)
    #print length
    #try 10 percent 
    percent_interval=(float(length)/100) *float(percent)
    #print int(percent_interval)
    
    #lowest=min(var_list)
    #highest=max(var_list)
    # print lowest, highest
    
    ## try to find intervals for truncation
    try_list=copy.deepcopy(var_list)
    try_list.sort()
    
    Thresholds = []
    # print list(chunks(try_list, int(percent_interval)))
    
    for x in list(chunks(try_list, int(percent_interval))):
     # print x[-1]
      Thresholds.append(x[-1])
    
    return  Thresholds

###End fly_threshold

################################
def Run_Rosetta(string, no_of_files, radius, Rosetta_cluster, RDB):  # rosetta can fail so repeat job
  rosettaout = os.getcwd()
  condition = 0

  while condition == 0:
    os.system (Rosetta_cluster+ ' -database '+RDB+' -cluster:radius ' +str(radius) + ' -in:file:s '+string+' > collect_out')

    my_sdtout = open (rosettaout + '/collect_out') 
    for line in my_sdtout:
     file_pattern = re.compile('^Clustering\s*(\d*)\s*structures')
     file_result = file_pattern.match(line)
     if file_result:
      #print line

      file_result2 = re.split(file_pattern, line )
      #print file_result2
      if int(file_result2[1]) !=  no_of_files:
       print 'fail there are ' + str(no_of_files) + ' files, got ' + file_result2[1]
      if int(file_result2[1]) ==  no_of_files:
       print ' SUCCESS there are ' + str(no_of_files) + ' files, got ' + file_result2[1] 
       condition = 1
####################################
def try_theseus(cmd):
    """
    #  theseus can fail so try Run a command with a timeout after which it will be forcibly killed.
    """
    has_worked = False
    while has_worked == False:
      Theseuslog=open(os.getcwd()+'/t.log', "w")
      p = subprocess.Popen(cmd  , shell =True, stdout =PIPE, stderr=PIPE)
    
    
      time.sleep(5)
      #print p.poll()
      if p.poll() is None: #still running      
         p.kill()
         #print 'timed out'
      else:
        #print p.communicate()
        has_worked = True
###END try_theseus

def Align_rosetta_fine_clusters_with_theseus(Rad_path, THESEUS, ):
    """
    If a cluster is present: align cluster  if no cluster is present: try to align most similar
    INPUT:
    Rad_path: directory with clustered PDB files, named C.0.X.pdb
    
    OUTPUT:
    Generates files named cluster_X..., but doesn't actually process them
    
    Theseus run with arguments:
     -a0  include alpha carbons and phosphorous atoms in superposition
     -r root name for output files
    """
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
    
    #print no_of_files
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

###END Align_rosetta_fine_clusters_with_theseus

def make_ensembles(trunc_out, threshold, THESEUS, MAX ):
    """
    Given a directory of truncated PDB files, use maxcluster to cluster them
    according to the three radius thresholds in RADS
    
    If there are more than 2 clusters returned by maxcluster copy them to the '/fine_clusters_'+str(RAD)
    directory, named as C.0.X, where X is the number of the cluster
    
    Theseus is then used to align the clusters and the path of the file containing the cluster is appended
    to the list ensembles_made, which is returned
    
    INPUTS:
    trunc_out: directory with truncated PDB files (e.g. fine_cluster_2/trunc_files_2)
    threshold: a float with the threshold the files were truncated with - used for file naming
    THESUS: path to theseus executable
    MAX: path to maxcluster executable
    
    """

    ensembles_made = []
    
    RADS = [1,2,3] # radius thresholds
    no_files = 0
    string = ''
    for infile in glob.glob( os.path.join(trunc_out, '*.pdb') ):
       string  = string + infile + ' '
       no_files +=1
    
    for RAD in RADS:
      Rad_path = trunc_out + '/fine_clusters_'+str(RAD)
      os.system('mkdir ' +trunc_out + '/fine_clusters_'+str(RAD))
      os.chdir(trunc_out + '/fine_clusters_'+str(RAD))
    
      #print 'no_files to sub cluster', no_files
      #print string, RAD, MAX, no_files
      cluster_files = cluster_with_MAX.cluster_with_MAX_FAST(string, RAD, MAX, no_files)  # use fastest method
    
      
      if cluster_files > 2:
       temp_name = 0
       for each_file in cluster_files:
         os.system('cp ' + each_file+' '+trunc_out + '/fine_clusters_'+str(RAD)+'/C.0.' + str(temp_name)+'.pdb' )
         temp_name+=1
      #Run_Rosetta(string, no_files,  RAD, Rosetta_cluster, RDB)
      
      ensemble_path = trunc_out + '/fine_clusters_'+str(RAD)+'_ensemble'
      os.system('mkdir '+ensemble_path )
      os.chdir(ensemble_path)
      
      # Run theseus to generate a  
      Align_rosetta_fine_clusters_with_theseus(Rad_path, THESEUS, )
      # Check if theseus worked and if so rename the file with the aligned files and append the path to the ensembles
      if os.path.exists(ensemble_path + '/cluster_0_sup.pdb'):
        os.system('mv ' +ensemble_path + '/cluster_0_sup.pdb '+ensemble_path + '/trunc_'+str(threshold)+'_rad_' +str(RAD)+'.pdb')
      ensembles_made.append(ensemble_path + '/trunc_'+str(threshold)+'_rad_' +str(RAD)+'.pdb')
    
    #print ensembles_made
    
    return ensembles_made
###END make_ensembles

def truncate(THESEUS, models_path, out_path, MAX, percent,FIXED_INTERVALS ):
    """
    Truncate the models in one folder
    * Run theseus to find the variances
    * For each truncation level, create a directory containing the truncated PDB files
    * Use maxcluster to cluster the truncated PDB
    
    INPUTS:
    THESUS: path to theseus
    models_path: RunDir + '/S_clusters/cluster_'+str(samples) [where clusters from spicker run are placed]
    out_path: RunDir+'/fine_cluster_'+str(samples) [ root directory where truncated clusters will go ]
    MAX: maxcluster executable
    percent: (var_args['percent'] or 5 ) - mutually exclusive with FIXED_INTERVALS. Truncate models so that "percent"
             residues are kept each cycle
    FIXED_INTERVALS: use the hard-coded intervals for truncation: [1, 1.5, 2 , 2.5, 3, 3.5 ,4, 4.5, 5, 5.5, 6, 7 , 8]
    
    OUTPUTS:
    truncate_log: file listing how many residues kept under each truncation threshold
    List of PDB files in trunc_out (e.g. fine_cluster_2/trunc_files_2.5) that contain the PDBs truncated to this level
    """

    all_ensembles = []
    
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
    #print run_dir
    
    os.system (THESEUS + ' -a0 `ls *.pdb | xargs`  >theseus_data')
    #print 'done'
    os.system('mv theseus* '+run_dir)
    #print run_dir
    T_data = out_path + '/theseus_variances.txt'
    
    #jmht - pause as it sometimes seem to take a while
    time.sleep(2)
    
    # for alternate versions of theseus remove RES cards
    # jmht - should probably just use a regular expression when we are processing the file
    tmp=open(out_path+'/tmp', "w")
    for line in open(out_path + '/theseus_variances.txt'):
        line = re.sub('RES ', '', line)
        tmp.write(line)
    tmp.close()
    os.system('mv '+out_path + '/theseus_variances.txt  '+  out_path + '/theseus_variances.txt_BAK')
    os.system('mv '+out_path+'/tmp '+  out_path + '/theseus_variances.txt')  
    
    #--------------------------------
    # choose threshold type
    #-------------------------------
    
    #FIXED_INTERVALS = False
    
    if FIXED_INTERVALS:
      thresholds = [1, 1.5, 2 , 2.5, 3, 3.5 ,4, 4.5, 5, 5.5, 6, 7 , 8]
    else:
      thresholds =fly_threshold(T_data, percent)
    
    # print thresholds, len(thresholds)
    
    
    #-------------------------------
    #truncate
    #----------------------------------
    truncate_log = open(out_path + '/truncate_log', "w")
    truncate_log.write('This is the number of residues kept under each truncation threshold\n\nthreshold\tnumber of residues\n')
    for threshold in thresholds:
      T_data = out_path + '/theseus_variances.txt'
      
      # Get a list of the indexes of the residues to keep
      add_list = make_list_to_keep(T_data, threshold)
      truncate_log.write( str(threshold) +'\t' + str(len(add_list)) + '\n' )
      trunc_out=''
    
      continue_trunc = True
      if len(add_list) < 1:  #######   limit to size of pdb to keep, if files are too variable, all residues are removed
         
         continue_trunc = False
    
    
      if continue_trunc == True:
       print 'truncating at '+str(threshold)
       trunc_out = out_path + '/trunc_files_' + str(threshold)
       #print trunc_out
       os.system('mkdir ' +trunc_out)
       
       for infile in glob.glob( os.path.join(models_path,  '*.pdb') ):
    
        name = re.split('/', infile)
        pdbname = str(name.pop())
    
        my_infile = open(infile, "r")
        pdbout = trunc_out+'/' + pdbname
        pdb_out = open (pdbout , "w")
        
        # Loop through PDB files and create new ones that only contain the residues left after truncation
        # Place the truncated PDB files in trunc_out
        for pdbline in my_infile:
          pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
          pdb_result = pdb_pattern.match(pdbline)
          if pdb_result:
            pdb_result2 = re.split(pdb_pattern, pdbline )
    
            for i in add_list : #convert to ints to compare
    
              if int(pdb_result2[5]) == int(i):
    
                 pdb_out.write(pdbline)
        
        pdb_out.close()
       #print 'making ensembles', threshold
       made_ens = make_ensembles(trunc_out, threshold, THESEUS, MAX )  ### MAKE ensemble for this trucnation level
       for a_ens in made_ens:
         all_ensembles.append(a_ens)
    #print 'got ensembles', all_ensembles 
    return all_ensembles

###END truncate

def truncate_Phenix(PHENIX, models_path, out_path, MAX, percent,FIXED_INTERVALS ):
  print 'asembling with Phenix'
  print PHENIX
  phenix_name = PHENIX

  if not re.search('ensembler', PHENIX):
      phenix_name = PHENIX+'.ensembler'
  print phenix_name

  #print 'here'
  all_ensembles = []

  run_dir = os.getcwd()
  if not os.path.exists(out_path):
   os.mkdir(out_path)
  os.chdir(out_path)
 
  number_of_models = 0
  list_of_pdbs = []
  string = ''

  os.chdir(models_path)
  #print run_dir

  for infile in glob.glob( os.path.join(models_path, '*.pdb') ):
     number_of_models +=1
     string  = string + infile + ' '
     list_of_pdbs.append(infile)

  cmd = phenix_name+' '+string
  #print cmd

  os.chdir(out_path)
  os.system(cmd)



  all_ensembles.append(os.path.join(out_path,'ensemble_merged.pdb'))
  #print all_ensembles
  return all_ensembles

###END truncate_Phenix
  
###################
#/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus 
#/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0 



if __name__ =="__main__":

 #            /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0
 THESEUS =  '/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus'
 models_path='/home/jaclyn/BEE/new_case/olga/models'
 out_path ='/home/jaclyn/DOMAINS/testing/truncted'
 MAX = '/home/jaclyn/programs/maxcluster/maxcluster'
 percent = 5

 #list_of_ensembles = truncate_Phenix(THESEUS, models_path, out_path, MAX, percent, True )

 rundir = '/home/jaclyn/BEE/ample_test/RO'
 list_of_ensembles = ['/home/jaclyn/BEE/ample_test/MODEL.pdb', '/home/jaclyn/BEE/ample_test/MODEL (copy).pdb']
 
 One_model_only(list_of_ensembles, rundir)









