#!/usr/bin/python

#
#  divide into parallel, select one, compare to rest- measure cluster repeat untill largest
#sort and rewrite
import re
import os, glob
import sys
import subprocess
import time
import pp
#import threading
from multiprocessing import Process, Queue, JoinableQueue


##############
def RunRosetta(cmd_string):       

   os.system(cmd_string)




def Run():
      cmd_string = 'ls '
      proc=2

      threads = []
      counter = 0
      i=0
      while i<proc: 


        thread = Process(target=RunRosetta, args=(cmd_string,))        
        thread.start() 
        threads.append(thread)

        i+=1
      print i  
      for thread in threads:

        thread.join()
        print thread.pid


###################################
#path = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/models' 
#path_to_LGA = '/home/jaclyn/Desktop/cmd/casp/LGA_package'
#outpath = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/run/'
#nproc = 1
#NMODELS = 1000#

#levels = [20,25,30,40,50,60]


#path = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/models' 


#outpath = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/' 
#nproc = 8 
#NMODELS = 10

#cluster_results = RUN_CLUSTER(levels, path, path_to_LGA, outpath, nproc, NMODELS)

#print 'GOT', cluster_results
Run()
