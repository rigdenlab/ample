#!/usr/bin/python2.6

#
#  Version 0.1  
#  no speed buffs, no extra functionality, no optimisation
#  Basic Pipeline
#  for multi core desktops, not clusters

# issues:
# 1) if this script is killed, the rosetta will continue
# 2) Theseus can hang, this is killed after a timeout 
# 3) is clustered by all atom RMSD to make fine clusters (option to cluster by CA only i included)
# 4) ASU content is the number of search models placed by MrBUMP. -- need to set this manually

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


from add_sidechains_SCWRL import *
import cluster_entropy, truncateedit, truncateedit_MAX
import SCWRL_edit, Final_display_results

import run_mr_bump_shelx_parallel
###############

      


#-------------------------------------------
#get command line options
#-------------------------------------------

parser = argparse.ArgumentParser(description='Structure solution by abinitio modeling', prefix_chars="-")

parser.add_argument('-ROSETTA', metavar='ROSETTA_path', type=str, nargs=1,
                   help='path for Rosetta AbinitioRelax')

parser.add_argument('-RDB', metavar='ROSETTA_database', type=str, nargs=1,
                   help='path for Rosetta database')

parser.add_argument('-fragsexe', metavar='path to make_fragments.pl', type=str, nargs=1,
                   help='location of make_fragments.pl')

parser.add_argument('-Rosetta_cluster', metavar='path to Rosettas cluster', type=str, nargs=1,
                   help='location of rosetta cluster')


parser.add_argument('-fasta', metavar='fasta_file', type=str, nargs=1,
                   help='protein fasta file. (required)')

parser.add_argument('-name', metavar='priotein name', type=str, nargs=1,
                   help='name of protien in the format ABCD_ ')


parser.add_argument('-NProc', metavar='NoProcessors', type=int, nargs=1,
                   help='number of processers (default 1)')


parser.add_argument('-RunDir', metavar='run_directory', type=str, nargs=1,
                   help='directory to put files (default current dir)')



parser.add_argument('-SCWRL', metavar='path to scwrl', type=str, nargs=1,
                   help='pathway to SCWRL exe')



parser.add_argument('-LGA', metavar='path_to_LGA dir', type=str, nargs=1,
                   help='pathway to LGA folder (not the exe) will use the \'lga\' executable')


parser.add_argument('-MAX', metavar='Maxcluster exe', type=str, nargs=1,
                   help='')


parser.add_argument('-THESEUS', metavar='Theseus exe (required)', type=str, nargs=1,
                   help='')

parser.add_argument('-MTZ', metavar='MTZ in', type=str, nargs=1,
                   help='')



#convert args to dictionary
args = parser.parse_args()
print args


var_args = vars(args)


#Required commands:
if var_args['ROSETTA'] is None:
   ROSETTA_PATH = 'no_path_given'
if not var_args['ROSETTA'] is None: 
   ROSETTA_PATH = var_args['ROSETTA'][0]

if var_args['RDB'] is None:
   ROSETTA_DB = 'no_path_given'
if not var_args['RDB'] is None: 
   ROSETTA_DB = var_args['RDB'][0]


if var_args['fragsexe'] is None:
  Make_fragents_exe = ' '
if not var_args['fragsexe'] is None:
  Make_fragents_exe = var_args['fragsexe'][0] 

if var_args['Rosetta_cluster'] is None:
   ROSETTA_cluster = 'no_path_given'
if not var_args['Rosetta_cluster'] is None: 
   ROSETTA_cluster = var_args['Rosetta_cluster'][0]


if var_args['fasta'] is None:
    FASTA = 'no_fasta_given'
if not var_args['fasta'] is None:
    FASTA = var_args['fasta'][0] 

if var_args['name'] is None:
   PDB_code = 'ABCD_'
if not var_args['name'] is None:
   PDB_code =var_args['name'][0]

if var_args['NProc'] is None:
   NProc = 1
if not var_args['NProc'] is None:
   NProc = var_args['NProc'][0]


if var_args['RunDir'] is None:
   RunDir = os.getcwd()
if not var_args['RunDir'] is None:
   RunDir = var_args['RunDir'][0]



if var_args['SCWRL'] is None:
  SCWRL = 'none'
if not var_args['SCWRL'] is None:
  SCWRL =var_args['SCWRL'][0] 


if var_args['LGA'] is None:
  LGA = 'none'
if not var_args['LGA'] is None:
  LGA =var_args['LGA'][0] 

if var_args['MAX'] is None:
  MAX = 'none'
if not var_args['MAX'] is None:
  MAX =var_args['MAX'][0]



if var_args['THESEUS'] is None:
   THESEUS = ''
if not var_args['THESEUS'] is None:
   THESEUS = var_args['THESEUS'][0]

if var_args['MTZ'] is None:
   MTZ = ''
if not var_args['MTZ'] is None:
   MTZ = var_args['MTZ'][0]

print  ROSETTA_PATH
print  FASTA
print PDB_code
print  NProc 
print  RunDir
print  SCWRL
print  LGA
print MAX
print THESEUS

#---------------------------------------
#check for errors
#---------------------------------------
#check if got all the programs

if not os.path.exists(ROSETTA_PATH):
    print 'You need to give the path for Rosetta abinito exe'
    sys.exit()




if not os.path.exists(FASTA):
    print 'You need to give the path for the fasta'
    sys.exit()

if not os.path.exists(RunDir):
    print 'You need to give a run directory'
    sys.exit()

if not os.path.exists(SCWRL):
    print 'You need to give the path for SCWRL'
    sys.exit()
if not os.path.exists(LGA):
    print 'You need to give the path for LGA'
    sys.exit()

if not os.path.exists(THESEUS):
    print 'You need to give the path for THESEUS'
    sys.exit()

if len(PDB_code)>5 or len(PDB_code)<5:
   print 'name is wrong format, use format ABCD_'
   sys.exit()

if PDB_code[4] != '_':
   print 'name is wrong format, use format ABCD_'
   sys.exit()  

if not os.path.exists(MTZ):
    print 'need mtz or no MR will be carried out'
    sys.exit()

#--------------------------------------------
#give the run file to check
#---------------------------------------------

Run_params = open(RunDir +'/Params_used', "w")
print var_args
for print_user in var_args:
 if var_args[print_user] is not None:
  print print_user +' : ' + str(var_args[print_user][0])

  Run_params.write(print_user +' : ' + str(var_args[print_user][0]) + '\n')
Run_params.close()




#-----------------------------------
#Do The Modelling
#-----------------------------------




time_start=time.time()
RUNNING=open( RunDir +'/ROSETTA.log', "w")


print PDB_code
os.chdir(RunDir)
  
  
#####make frags 

RUNNING.write('----- making fragments--------\n')
RUNNING.flush()
frags_dir = RunDir + '/frags'
os.system('mkdir ' + frags_dir)
print Make_fragents_exe + ' -rundir ' + frags_dir + ' -id ' + PDB_code + ' '+FASTA+ ' -nojufo -nosam -noprof -nohoms ' 
os.system(Make_fragents_exe + ' -rundir ' + frags_dir + ' -id ' + PDB_code + ' '+FASTA+ ' -nojufo -nosam -noprof -nohoms '   ) 

frags_3_mers = frags_dir + '/aa' +PDB_code+'03_05.200_v1_3'
frags_9_mers = frags_dir + '/aa' +PDB_code+'09_05.200_v1_3'

RUNNING.write('Fragments done\n3mers at: '+frags_3_mers+'\n9mers at: '+frags_9_mers+'\n\n')

###model
RUNNING.write('----- making models--------\n')
RUNNING.flush()
DEBUG = False
NMODELS = 1000
MakeModels = True 
if MakeModels == True:
      
      previous_seeds = [0]    #make random seeds  (1000000, 6000000) ! Must be unique seeds!!!!
      proc = 1
      while proc < NProc +1:
        new_seed = random.randint(1000000, 6000000)  
        seed_present = False

        for seed in previous_seeds:
           if new_seed == seed:
             seed_present = True
             break

        if seed_present == False:
          previous_seeds.append(new_seed)
          proc +=1

      print previous_seeds


      split_jobs =  NMODELS / NProc   ### split jobs between processors
      remainder =   NMODELS % NProc
      jobs = [0]      
      proc = 1
      while proc < NProc +1:
        jobs.insert(proc, split_jobs)
        proc +=1
      jobs[-1] = jobs[-1] + remainder    
      print jobs  ##################################




      os.system('mkdir '+ RunDir +'/models')
      PATH_TO_MODELS = RunDir +'/models'

      proc = 1
      while proc < NProc +1:
       print proc
 
       os.system('mkdir '+ RunDir + '/models_'+str(proc))
       os.chdir(RunDir + '/models_'+str(proc))

       print  ROSETTA_PATH +' -database ' + ROSETTA_DB + ' -in::file::fasta ' + FASTA + ' -in:file:frag3 '+ frags_3_mers +' -in:file:frag9 '+ frags_9_mers + ' -out:path ' +RunDir + '/models_'+str(proc)+' -out:pdb -out:nstruct ' + str(jobs[proc]) + ' -out:file:silent '+RunDir +'models_1/OUT -return_full_atom false   -run:constant_seed -run:jran ' + str( previous_seeds[proc]) + '  &'

       RUNNING.write(ROSETTA_PATH +' -database ' + ROSETTA_DB + ' -in::file::fasta ' + FASTA + ' -in:file:frag3 '+ frags_3_mers +' -in:file:frag9 '+ frags_9_mers + ' -out:path ' +RunDir + '/models_'+str(proc)+' -out:pdb -out:nstruct ' + str(jobs[proc]) + ' -out:file:silent '+RunDir +'models_1/OUT -return_full_atom false   -run:constant_seed -run:jran ' + str( previous_seeds[proc]) + '  &\n')
       RUNNING.flush()


       os.system(ROSETTA_PATH +' -database ' + ROSETTA_DB + ' -in::file::fasta ' + FASTA + ' -in:file:frag3 '+ frags_3_mers +' -in:file:frag9 '+ frags_9_mers + ' -out:path ' +RunDir + '/models_'+str(proc)+' -out:pdb -out:nstruct ' + str(jobs[proc]) + ' -out:file:silent '+RunDir +'/models_'+str(proc) +'/OUT -return_full_atom false   -run:constant_seed -run:jran ' + str( previous_seeds[proc]) + '  &')

       proc +=1




### wait for all models to be made, check number for each job

      no_models_to_make = NMODELS
      no_models_have = 0 
      while no_models_have !=  no_models_to_make:
        no_models_have = 0
        proc = 1
        while proc < NProc +1:
          list_of_files = [file for file in os.listdir(RunDir +'/models_'+str(proc) ) if file.lower().endswith('.pdb')]


          no_models_have  += len(list_of_files)
          proc+=1
        print 'Number of models made so far = ' + str(no_models_have)
        if no_models_have ==  no_models_to_make:
           break
        time.sleep(5)       
   
     
   

      
   
      LOG = open(RunDir+'/LOG', "a")
      print 'MODELLING DONE'
## got them all   Must Delete Models or job will hang -----------



      MODELLING_time_stop=time.time()
      cat_string = ''

     # RunDir +'/models'
      proc = 1
      while proc < NProc +1:
         print proc
         
         add_sidechains_SCWRL(SCWRL,   RunDir + '/models_'+str(proc),   RunDir + '/models',  str(proc) )
         cat_string += RunDir + '/OUT_' + str(proc) +'  '
         #os.system('rm -r '+RunDir + '/models_'+str(proc) )
         proc+=1
    
      os.chdir(RunDir)
      print '\n' +cat_string


     ####### keep or delete run files


PATH_TO_MODELS = RunDir + '/models'
RUNNING.write('\n\nModelling done\nmodels at '+PATH_TO_MODELS+'\n\n')
#--------------------------------------
# Do the clustering  
#---------------------------------------

levels = [20,25,30,40,50,60]   # clustering levels
print '----------CLUSTERING---------'
RUNNING.write('----------CLUSTERING---------\n')
RUNNING.write('Default clustering to find only the largest cluster (see developemnt version for alternate clustering methods)')
RUNNING.flush()
cluster_path = RunDir + '/clusters/'  # outpath for clusters
os.system('mkdir '+ RunDir + '/clusters')



cluster_results = cluster_entropy.RUN_CLUSTER(levels, PATH_TO_MODELS, LGA, cluster_path,  NProc, NMODELS)


print 'Clustering Done!'
RUNNING.write('Clustering Done! your results:\n')
i = 0
while i < len(cluster_results):
 RUNNING.write('at level ' +str(levels[i]) + ' cluster size is ' + str(cluster_results[i])+'\n' )
 i+=1
RUNNING.flush()
#----------------------------------------
#pick a clustering level
#----------------------------------------
#  default method of picking a cluster  cutoffs <50 >5:



trunc_level = 'nan'
print levels, cluster_results
i = 0
while i < len(cluster_results):
 if  cluster_results[i] <50 and cluster_results[i] >5:
  trunc_level = levels[i]
  break
 i+=1
print trunc_level


#if none exist in range, go higher:
if trunc_level == 'nan':
  i = 0
  while i < len(cluster_results):
    if  cluster_results[i] >5 :
      trunc_level = levels[i]
    i+=1

print trunc_level

#-------------------------------------
#Truncate
#-------------------------------------
RUNNING.write('\n----------Truncating---------\n')
RUNNING.flush()
os.system('mkdir '+RunDir+'/fine')
os.chdir( RunDir+'/fine')
models_path = RunDir + '/clusters/cluster_'+str(trunc_level)+'/sorted_cluster_0'

print THESEUS, models_path, RunDir+'/fine', ROSETTA_cluster 

Strict= True #  Default = True
if Strict== True:
 list_of_ensembles = truncateedit.truncate(THESEUS, models_path, RunDir+'/fine', ROSETTA_cluster, ROSETTA_DB )
if Strict== False:
  list_of_ensembles = truncateedit_MAX.truncate(THESEUS, models_path, RunDir+'/fine', MAX )

RUNNING.write('Truncating done!\n')
RUNNING.flush()
#-------------------------------------
#fix sidechains
#------------------------------------
os.system('mkdir '+RunDir+'/ensembles')


for each_ens in list_of_ensembles:
 SCWRL_edit.edit_sidechains(each_ens, RunDir+'/ensembles/')


final_ensembles =[]
for infile in glob.glob( os.path.join(RunDir+'/ensembles', '*.pdb') ):
  final_ensembles.append(infile)

#------------------------------------
# END
#----------------------------------
time_stop=time.time()

elapsed_time= time_stop - time_start
run_in_min=elapsed_time/60
run_in_hours = run_in_min/60

RUNNING.write('\nMODELLING and ASSEMBLY ALL DONE  (in '+str(run_in_hours)+' hours) \n----------------------------------------\n')
RUNNING.flush()



#--------------Import into MrBUMP here-------------------
RUNNING.write('Running MrBUMP\nMR LOGS in '+RunDir+'/MRBUMP')
RUNNING.flush()

os.system('mkdir ' + RunDir+'/MRBUMP')
os.chdir(RunDir+'/MRBUMP')
bump_dir = RunDir+'/MRBUMP'

split_ensembles = run_mr_bump_shelx_parallel.split(final_ensembles, NProc)
run_mr_bump_shelx_parallel.split_into_runs(MTZ, split_ensembles, bump_dir, FASTA, NProc)
Final_display_results.make_log(os.path.join(bump_dir, RunDir, 'Final_results.log'))

time_stop=time.time()
elapsed_time= time_stop - time_start
run_in_min=elapsed_time/60
run_in_hours = run_in_min/60

RUNNING.write('\nMR and shelx ALL DONE  (in '+str(run_in_hours)+' hours) \n----------------------------------------\n')
RUNNING.flush()

