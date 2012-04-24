#!/usr/bin/python2.6

#edit the sidechains to make polyala, all and reliable

import re
import os, glob
import sys

from itertools import repeat

###########################################
def get_clusters_from_distances(matrix, radius):

  #print len(matrix)
  matrix_line = 0
  best_cluster=[]
  largest = 0
 

  while matrix_line<len(matrix):
   cluster=0
   current_cluster_models = []
   
   each_model = 0

   while each_model < len(matrix[matrix_line]):
     if float(matrix[matrix_line][each_model])<radius:
       current_cluster_models.append(each_model)
       cluster+=1
     each_model+=1
  # print 'MODEL ', matrix_line+1, cluster
   if len(current_cluster_models) > largest:
      largest = len(current_cluster_models)
      best_cluster = current_cluster_models

   matrix_line +=1
   
 # print 'BEST ', best_cluster
  return best_cluster
########################################## ADD ALL
def cluster_with_MAX_FASTBAK(string, radius, MAX, no_models):
  cur_dir = os.getcwd()
  string = re.sub(' ', '\n', string)
  list_string = open(cur_dir + '/list', "w")
  list_string.write(string)
  list_string.close()  

  model_indeces=[]
  #print 'runing MAX'
  os.system(MAX + ' -l list  -L 4 -rmsd -d 1000 -bb -C0 >MAX_LOG ')
  #print 'MAX Done'  
  matrix = []
  inc = 1
  while inc < no_models+1:
      matrix_line =[]
      #print 'matrix', inc
      max_log = open(cur_dir+'/MAX_LOG')
      pattern = re.compile('INFO  \: Model')
      for line in max_log:
        if re.match(pattern, line):

         #print line 
         split = re.split('INFO  \: Model\s*(\d*)\s*(.*)\.pdb\s*vs\. Model\s*(\d*)\s*(.*)\.pdb\s*=\s*(\d*\.\d*)', line)
         #print split 
         if int(split[1]) == inc:
           matrix_line.insert(int(split[3]),  split[5] )
           if split[2] not in model_indeces: 
            model_indeces.insert(inc, split[2])

         if int(split[3]) == inc:
           matrix_line.insert(int(split[1]), split[5] )
           if split[4] not in model_indeces: 
            model_indeces.insert(inc, split[4])


      matrix_line.insert(inc-1, 0 )
      matrix.append(matrix_line ) 
      inc +=1

 # for x in matrix:   ### got distance matrix
 #  print x
 # print 'GOT MATRIX'
  cluster_models = get_clusters_from_distances(matrix, radius)
  return cluster_models
#####################################
##################################### START

def matrix_insert(i, j, k,  matrix):
 matrix[1][2]  = 2
 #print matrix
 
 return matrix

##########
def cluster_with_MAX_FAST(string, radius, MAX, no_models):
  cur_dir = os.getcwd()
  string = re.sub(' ', '\n', string)
  list_string = open(cur_dir + '/list', "w")
  list_string.write(string)
  list_string.close()  

  models=[0]*no_models
  #print 'runing MAX'
  os.system(MAX + ' -l list  -L 4 -rmsd -d 1000 -bb -C0 >MAX_LOG ')
  #print 'MAX Done'  


  matrix =  [[0 for col in range(no_models)] for row in range(no_models)]
  max_log = open(cur_dir+'/MAX_LOG')
  pattern = re.compile('INFO  \: Model')
  for line in max_log:
        if re.match(pattern, line):


         split = re.split('INFO  \: Model\s*(\d*)\s*(.*)\.pdb\s*vs\. Model\s*(\d*)\s*(.*)\.pdb\s*=\s*(\d*\.\d*)', line)
         #print split 
#         int(split[3])

        # print split[1], split[3], split[5]
         matrix[  int(split[1]) -1 ][  int(split[3]) -1]  = split[5]

         if split[2]+'.pdb' not  in models:
                models[int(split[1]) -1]  =  split[2]+'.pdb'

         if split[4]+'.pdb' not  in models:
                models[int(split[3]) -1]  =  split[4]+'.pdb'

  x = 0
  while x < len(matrix):
    y = 0
    while y < len(matrix):
     matrix[y][x] = matrix [x][y]
     y+=1
    x+=1

 # for x in matrix:   ### got distance matrix
 #  print x
 # print 'GOT MATRIX'
############# get model indeces


  CLUSTER = []

  cluster_models = get_clusters_from_distances(matrix, radius)
  #print cluster_models
  #for i in cluster_models:
  #  print 'using', i

  for index in cluster_models:
    CLUSTER.append(models[index])
  #print CLUSTER
  return CLUSTER




#  1 2 3 4
#1
#2
#3
#4




#################### STOP
def cluster_with_MAX(string, radius, MAX, no_models):
  cur_dir = os.getcwd()
  string = re.sub(' ', '\n', string)
  list_string = open(cur_dir + '/list', "w")
  list_string.write(string)
  list_string.close()  



  os.system(MAX + ' -l list  -L 4 -rmsd -d 1000 -bb >MAX_LOG ')
  
  matrix = []
  inc = 1

  largest_size = 0
  largest_models = []

  while inc < no_models+1:  # for each model
      COUNT = 1   # size of cluster
      models = []

      max_log = open(cur_dir+'/MAX_LOG')
      pattern = re.compile('INFO  \: Model')
      for line in max_log:
        if re.match(pattern, line):
       #  print line 
         split = re.split('INFO  \: Model\s*(\d*)\s*(.*)\.pdb\s*vs\. Model\s*(\d*)\s*(.*)\.pdb\s*=\s*(\d*\.\d*)', line)
       #  print split




         if int(split[1]) == inc:
        
           if float(split[5]) < radius:
              COUNT +=1
              models.append(split[4]+'.pdb')
              if split[2] not  in models:
                models.append(split[2]+'.pdb')              

         if int(split[3]) == inc:
           if float(split[5]) < radius:
              COUNT +=1
              models.append(split[2]+'.pdb')
              if split[4] not  in models:
                models.append(split[4]+'.pdb')  


      if COUNT > largest_size:
         largest_size = COUNT
         largest_models = models
     # print 'MODEL ',  inc,  ' ', COUNT
    #  for x in models:
     #     print x

      inc +=1

 # print 'LARGEST ',  largest_size
  #for l in largest_models:
  #        print l
  return largest_models

###################
#string ='/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/8_S_00000002.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/8_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/2_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/5_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/1_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/3_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/8_S_00000003.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/4_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/7_S_00000001.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/fine/trunc_files_1/6_S_00000001.pdb' 
#1 /home/jaclyn/programs/maxcluster/maxcluster 10

#string = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000009.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000010.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000011.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000012.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000013.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000016.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000017.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000018.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000019.pdb /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0/1_S_00000020.pdb '
#NoMODELS = 10
#
#cluster_with_MAX_FAST(string, 2, '/home/jaclyn/programs/maxcluster/maxcluster', NoMODELS)

