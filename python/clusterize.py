#!/usr/bin/env python

# Class for submitting ample modelling and MR/Refine/Build components to a cluster
# queuing system.
#
# Ronan Keegan 25/10/2011
#

import os
import re
import sys
import subprocess
import shlex
import string
import time
import random
import run_shelx
import MTZParse
import MatthewsCoef

if not "CCP4" in sorted(os.environ.keys()):
    raise RuntimeError('CCP4 not found')

###############

class ClusterRun:

   def __init__(self):
   
      self.qList=[]
      self.runningQueueList=[]
      self.QTYPE=""
      self.SCWRL_EXE=""
      self.USE_SCWRL=False

      self.RunDir=""
      self.jobLogsList=[]
   
      self.SEQIN=""
      self.HKLIN=""
      self.LABIN=dict([])
      self.MTZcell=[]
      self.MTZspacegroup=""
      self.MTZresolution=0.0
      self.MTZcolLabels=dict([])
      
      self.BCYCLES  = 5
      self.BUCC = True
      self.ACYCLES  = 5
      self.ARPWARP = True
      self.SHELXE = False
      self.SCYCLES = 15


      #self.shelxClusterScript="python " + os.path.join(os.environ["CCP4"], "share", "ample", "python", "shelx_cluster.py")
      self.shelxClusterScript="python " + os.path.join(os.environ["CCP4"], "share", "ample", "python", "shelxe_trace.py")

      self.ALLATOM=False
      self.debug=True

      if os.name == "nt":
         self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset.exe")
      else:
         self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset")

   def set_USE_SCWRL(self, bool):
      self.USE_SCWRL=bool
      
   def setScwrlEXE(self, exePath):
      self.SCWRL_EXE=exePath

   def getMTZInfo(self, HKLIN, workingDir):
      """ Extract useful information from an MTZ file """

      mtz_info=MTZParse.Mtzdump()

      mtz_info.setHKLIN(HKLIN)
      mtz_info.setMTZdumpLogfile(workingDir) 
      mtz_info.go() 
      mtz_info.getColumnData()

      self.MTZcell=mtz_info.getCell()
      self.MTZspacegroup=mtz_info.getSpacegroup()
      if string.strip(self.MTZspacegroup.replace(" ", ""))=="P-1":
         self.MTZspacegroup="P 1"
      self.MTZresolution=float(mtz_info.getResolution())
      self.MTZcolLabels=mtz_info.colLabels
       


   def setupModellingDir(self, RunDir):
      """ A function to create the necessary directories for the modelling step """
  
      self.RunDir=RunDir

      # Create the directories for the submission scripts
      if os.path.isdir(os.path.join(RunDir, "models")) == False:
         os.mkdir(os.path.join(RunDir, "models"))
      if os.path.isdir(os.path.join(RunDir, "pre_models")) == False:
         os.mkdir(os.path.join(RunDir, "pre_models"))
      if os.path.isdir(os.path.join(RunDir, "pre_models", "sge_scripts")) == False:
         os.mkdir(os.path.join(RunDir, "pre_models", "sge_scripts"))
      if os.path.isdir(os.path.join(RunDir, "pre_models", "logs")) == False:
         os.mkdir(os.path.join(RunDir, "pre_models", "logs"))


   def monitorQueue(self, user=""):
      """ Monitor the Cluster queue to see when all jobs are completed """

      sys.stdout.write("Jobs submitted to cluster queue, awaiting their completion...\n") 

      # set a holder for the qlist 
      runningList=self.qList
      newRunningList=[]

      while runningList!=[]:
         time.sleep(60)
         self.getRunningJobList(user)
         for job in runningList: 
            if str(job) in self.runningQueueList:
               newRunningList.append(job)
         if len(runningList) > len(newRunningList):
            sys.stdout.write("Queue Monitor: %d out of %d jobs remaining in cluster queue...\n" %  (len(newRunningList),len(self.qList)))   
         if len(newRunningList) == 0:
            sys.stdout.write("Queue Monitor: All jobs complete!\n")   
         runningList=newRunningList
         newRunningList=[]
         
   def getJobStatus(self, qNumber):
      """ Check a job status int the cluster queue """
 
      status=1 

      command_line='qstat -j %d' % qNumber

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.PIPE)
  
      (child_stdout, child_stderr, child_stdin) = (p.stdout, p.stderr, p.stdin)
  
      # Write the keyword input
      child_stdin.close()
  
      child_stdout.close()

      # Watch the output for successful termination
      err=child_stderr.readline()
  
      while err:
         #sys.stdout.write(out)
         if self.QTYPE=="SGE":
            if "Following jobs do not exist" in err:
               status=0
         err=child_stderr.readline()
  
      child_stderr.close()

      return status

   def getRunningJobList(self, user=""):
      """ Check a job status int the cluster queue """
 
      if user == "":
         command_line='qstat'
      else:
         command_line='qstat -u ' + user

      log_lines=[]
      self.runningQueueList=[]

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE)
  
      (child_stdout, child_stdin) = (p.stdout, p.stdin)
  
      # Write the keyword input
      child_stdin.close()
  
      # Read the output
      out=child_stdout.readline()
      
      while out:
         #sys.stdout.write(out)
         if self.QTYPE=="SGE":
            log_lines.append(string.strip(out))
         out=child_stdout.readline()
  
      child_stdout.close()
   
      if log_lines != []:
         log_lines.pop(0)
         log_lines.pop(0)
         for i in log_lines:
            self.runningQueueList.append(string.split(i)[0])

   def NMRmodelOnCluster(self, RunDir, proc, jobNumber, ROSETTA_PATH, ROSETTA_DB, FASTA, frags_3_mers, frags_9_mers, ideal_homolog,  ALI, seed, MR_ROSETTA ):
      """ Farm out the modelling step on a cluster (SGE) """
      # Set the file number according to the job number
      if jobNumber<10:
         fileNumber="0000000" + str(jobNumber)
      elif jobNumber>=10 and jobNumber<100:
         fileNumber="000000" + str(jobNumber)
      elif jobNumber>=100 and jobNumber<1000:
         fileNumber="00000" + str(jobNumber)
      elif jobNumber>=1000 and jobNumber<10000:
         fileNumber="0000" + str(jobNumber)
      elif jobNumber>=10000 and jobNumber<100000:
         fileNumber="000" + str(jobNumber)
      elif jobNumber>=100000 and jobNumber<1000000:
         fileNumber="00" + str(jobNumber)
      else:
         sys.stdout.write("No. of Models exceeds program limits (Max=999999)\n")
         sys.exit()
      preModelDir=os.path.join(RunDir, "pre_models", "model_" + str(jobNumber))

      if os.path.isdir(preModelDir) == False:
         os.mkdir(preModelDir)

      n = os.path.split(ideal_homolog)
      n ='S_' +n[-1].rstrip('.pdb')+'_0001.pdb'
        

      PDBInFile     = os.path.join(preModelDir, n)
      PDBSetOutFile = os.path.join(preModelDir, "pdbsetOut_" + str(jobNumber) + ".pdb")
      PDBScwrlFile  = os.path.join(preModelDir, "scwrlOut_" + str(jobNumber) + ".pdb")
      SEQFile       = os.path.join(preModelDir, "S_" + fileNumber + ".seq")
      PDBOutFile    = os.path.join(RunDir, "models", "1_S_" + fileNumber + ".pdb")

     # Create a cluster submission script for this modelling job
      jobName="model_" + str(proc) + "_" + str(seed)
      sub_script=os.path.join(RunDir, "pre_models", "sge_scripts", "job_" + jobName + ".sub")

      self.jobLogsList.append(os.path.join(RunDir, "pre_models", "logs", jobName + '.log'))

      file=open(sub_script, "w")
      file.write('#!/bin/sh\n'
         '#$ -j y\n' +
         '#$ -cwd\n' +
         '#$ -w e\n' +
         '#$ -V\n' +
         '#$ -o ' + os.path.join(RunDir, "pre_models", "logs", jobName + '.log') + '\n' +
         '#$ -N ' + jobName + '\n\n')

      file.write("setenv CCP4_SCR $TMPDIR\n\n")
      
      file.write('cd '+ os.path.join(RunDir, "pre_models", "model_" + str(jobNumber)) +'\n\n'+
           MR_ROSETTA +' \\\n'+
           '-database '+ROSETTA_DB+' \\\n'+
           '-MR:mode cm \\\n'+
           '-in:file:extended_pose 1 \\\n'+
           '-in:file:fasta '+FASTA+' \\\n'+
           '-in:file:alignment '+ALI+' \\\n'+
           '-in:file:template_pdb '+ideal_homolog+' \\\n'+
           '-loops:frag_sizes 9 3 1 \\\n'+
           '-loops:frag_files '+frags_9_mers+' '+frags_3_mers+' none \\\n'+
           '-loops:random_order \\\n'+
           '-loops:random_grow_loops_by 5 \\\n'+
           '-loops:extended \\\n'+
           '-loops:remodel quick_ccd \\\n'+
           '-loops:relax relax \\\n'+
           '-relax:default_repeats 4 \\\n'+
           '-relax:jump_move true    \\\n'+
           '-cm:aln_format grishin \\\n'+
           '-MR:max_gaplength_to_model 8 \\\n'+
           '-nstruct 1  \\\n'+
           '-ignore_unrecognized_res \\\n'+
           '-overwrite \n\n')
           

      file.write("pushd " + os.path.join(preModelDir) + "\n\n" +

      self.pdbsetEXE + " xyzin " + PDBInFile + " xyzout " + PDBSetOutFile + "<<eof\n" +
      "sequence single\n" +
      "eof\n\n" +

      "tail -n +2 SEQUENCE | sed s'/ //g' >> " + SEQFile + "\n" +
      "popd\n\n"  )
      if self.USE_SCWRL :

          file.write(self.SCWRL_EXE + " -i " + PDBInFile + " -o " + PDBScwrlFile + " -s " + SEQFile + "\n\n" +
          "head -n -1 " + PDBScwrlFile + " >> " + PDBOutFile + "\n" +
           "\n")
      if not self.USE_SCWRL :
          file.write('cp ' + PDBInFile + ' ' +  PDBOutFile + "\n" )

      # Clean up non-essential files unless we are debugging
      if self.debug == False:
         file.write("rm " + PDBSetOutFile + "\n" +
         "rm " + os.path.join(preModelDir, "SEQUENCE") + "\n" +
         "rm " + PDBScwrlFile + "\n\n")

      file.close()
      
      # Submit the job
      curDir=os.getcwd()
      os.chdir(os.path.join(RunDir, "pre_models", "sge_scripts"))
      command_line='qsub -V %s' % sub_script

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                    stdout = subprocess.PIPE)

      (child_stdout, child_stdin) = (p.stdout, p.stdin)

      # Write the keyword input
      child_stdin.close()

      # Watch the output for successful termination
      out=child_stdout.readline()

      qNumber=0

      while out:
         #sys.stdout.write(out)
         if self.QTYPE=="SGE":
            if "Your job" in out:
               qNumber=int(string.split(out)[2])
               self.qList.append(qNumber)
         out=child_stdout.readline()

      child_stdout.close()

      os.chdir(curDir)



   def modelOnCluster(self, RunDir, proc, jobNumber, ROSETTA_PATH, ROSETTA_DB, FASTA, frags_3_mers, frags_9_mers, seed, rosetta_string):
      """ Farm out the modelling step on a cluster (SGE) """

      # Set the file number according to the job number
      if jobNumber<10:
         fileNumber="0000000" + str(jobNumber)
      elif jobNumber>=10 and jobNumber<100:
         fileNumber="000000" + str(jobNumber)
      elif jobNumber>=100 and jobNumber<1000:
         fileNumber="00000" + str(jobNumber)
      elif jobNumber>=1000 and jobNumber<10000:
         fileNumber="0000" + str(jobNumber)
      elif jobNumber>=10000 and jobNumber<100000:
         fileNumber="000" + str(jobNumber)
      elif jobNumber>=100000 and jobNumber<1000000:
         fileNumber="00" + str(jobNumber)
      else:
         sys.stdout.write("No. of Models exceeds program limits (Max=999999)\n")
         sys.exit()

      preModelDir=os.path.join(RunDir, "pre_models", "model_" + str(jobNumber))
      
      if os.path.isdir(preModelDir) == False:
         os.mkdir(preModelDir)

      PDBInFile     = os.path.join(preModelDir, "S_00000001.pdb")
      PDBSetOutFile = os.path.join(preModelDir, "pdbsetOut_" + str(jobNumber) + ".pdb")
      PDBScwrlFile  = os.path.join(preModelDir, "scwrlOut_" + str(jobNumber) + ".pdb")
      SEQFile       = os.path.join(preModelDir, "S_" + fileNumber + ".seq")
      PDBOutFile    = os.path.join(RunDir, "models", "1_S_" + fileNumber + ".pdb")

      # Create a cluster submission script for this modelling job
      jobName="model_" + str(proc) + "_" + str(seed)
      sub_script=os.path.join(RunDir, "pre_models", "sge_scripts", "job_" + jobName + ".sub")

      self.jobLogsList.append(os.path.join(RunDir, "pre_models", "logs", jobName + '.log'))

      file=open(sub_script, "w")
      file.write('#!/bin/sh\n'
         '#$ -j y\n' +
         '#$ -cwd\n' +
         '#$ -w e\n' +
         '#$ -V\n' +
         '#$ -o ' + os.path.join(RunDir, "pre_models", "logs", jobName + '.log') + '\n' +
         '#$ -N ' + jobName + '\n\n')
  
      file.write("setenv CCP4_SCR $TMPDIR\n\n")
      
      if self.ALLATOM:
         file.write(ROSETTA_PATH +' -database ' + ROSETTA_DB + ' -in::file::fasta ' + FASTA +
         ' -in:file:frag3 '+ frags_3_mers +' -in:file:frag9 '+ frags_9_mers +
         ' -out:path ' + preModelDir +' -out:pdb -out:nstruct 1' +
         ' -out:file:silent '+ preModelDir +
         '/OUT -return_full_atom true -abinitio:relax -run:constant_seed -run:jran ' + str(seed) + rosetta_string + "\n\n")
      else:
         file.write(ROSETTA_PATH +' -database ' + ROSETTA_DB + ' -in::file::fasta ' + FASTA +
         ' -in:file:frag3 '+ frags_3_mers +' -in:file:frag9 '+ frags_9_mers +
         ' -out:path ' + preModelDir +' -out:pdb -out:nstruct 1' +
         ' -out:file:silent '+ preModelDir +
         '/OUT -return_full_atom false -run:constant_seed -run:jran ' + str(seed) + rosetta_string + "\n\n") 

      file.write("pushd " + os.path.join(preModelDir) + "\n\n" +
 
      self.pdbsetEXE + " xyzin " + PDBInFile + " xyzout " + PDBSetOutFile + "<<eof\n" +
      "sequence single\n" +
      "eof\n\n" +

      "tail -n +2 SEQUENCE | sed s'/ //g' >> " + SEQFile + "\n" + 
      "popd\n\n"  )
      if self.USE_SCWRL :
      
          file.write(self.SCWRL_EXE + " -i " + PDBInFile + " -o " + PDBScwrlFile + " -s " + SEQFile + "\n\n" +
          "head -n -1 " + PDBScwrlFile + " >> " + PDBOutFile + "\n" + 
           "\n")
      if not self.USE_SCWRL :
          file.write('cp ' + PDBInFile + ' ' +  PDBOutFile + "\n" )  
  
      # Clean up non-essential files unless we are debugging
      if self.debug == False:
         file.write("rm " + PDBSetOutFile + "\n" + 
         "rm " + os.path.join(preModelDir, "SEQUENCE") + "\n" +
         "rm " + PDBScwrlFile + "\n\n")

      file.close()
  
      # Submit the job
      curDir=os.getcwd()
      os.chdir(os.path.join(RunDir, "pre_models", "sge_scripts"))
      command_line='qsub -V %s' % sub_script
  
      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                    stdout = subprocess.PIPE)
  
      (child_stdout, child_stdin) = (p.stdout, p.stdin)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout.readline()
  
      qNumber=0
  
      while out:
         #sys.stdout.write(out)
         if self.QTYPE=="SGE":
            if "Your job" in out:
               qNumber=int(string.split(out)[2])
               self.qList.append(qNumber)
         out=child_stdout.readline()
  
      child_stdout.close()
  
      os.chdir(curDir)
  
   def mrBuildOnCluster(self, clusterDir, ensemblePDB, jobID,mrbump_programs, fixedPDB="", fixedIDEN="", ):
      """ Run the molecular replacement and model building on a cluster node """

      # Create a cluster submission script for this modelling job
      jobName="mrBuild_" + str(jobID) 

      jobDir=os.path.join(clusterDir, "mrbuild_" + str(jobID))
      os.mkdir(jobDir)
      os.mkdir(os.path.join(jobDir, "sge_scripts"))
      os.mkdir(os.path.join(jobDir, "logs"))
      sub_script=os.path.join(jobDir, "sge_scripts", "job_" + str(jobID) + ".sub")

      # Get details about mol weight and solvent content and number of molecules in the ASU
      matthews_coef=MatthewsCoef.MattCoef()

      CELLstring=" ".join(map(str, self.MTZcell))

      matthews_coef.setCELL(CELLstring)
      matthews_coef.setSYMM(self.MTZspacegroup)
      matthews_coef.setRESO(self.MTZresolution)

      MOLWT = float(run_shelx.seqwt(self.SEQIN))
      matthews_coef.runMC(MOLWT, os.path.join(jobDir, "logs", "matthes_coef.log"))

      if self.MTZresolution<2.0:
         free_lunch = ' -e1.0 -l8'
        
      else:
         free_lunch = ' '
        
      modelName=os.path.split(ensemblePDB)[1].replace(".pdb","")

      phaserPDB = os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data', 'loc0_ALL_'+ modelName +'', 'unmod', 'mr', 'phaser', 'refine', 'refmac_phaser_loc0_ALL_'+ modelName +'_UNMOD.pdb')
      phaserMTZ = os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data', 'loc0_ALL_'+ modelName +'', 'unmod', 'mr', 'phaser', 'refine', 'refmac_phaser_HKLOUT_loc0_ALL_'+ modelName +'_UNMOD.mtz')
      molrepPDB = os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data', 'loc0_ALL_'+ modelName +'', 'unmod', 'mr', 'molrep', 'refine', 'refmac_molrep_loc0_ALL_'+ modelName +'_UNMOD.pdb')
      molrepMTZ = os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data', 'loc0_ALL_'+ modelName +'', 'unmod', 'mr', 'molrep', 'refine', 'refmac_molrep_HKLOUT_loc0_ALL_'+ modelName +'_UNMOD.mtz')
       
      phaserTFZPDB = os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data', 'loc0_ALL_'+ modelName +'', 'unmod', 'mr', 'phaser', 'refine', 'phaser_loc0_ALL_'+ modelName +'_UNMOD.1.pdb')
     
      file=open(sub_script, "w")
      file.write('#!/bin/sh\n'
         '#$ -j y\n' +
         '#$ -cwd\n' +
         '#$ -w e\n' +
         '#$ -V\n' +
         '#$ -o ' + os.path.join(jobDir, "logs", "job_" + str(jobID) + '.log') + '\n' +
         '#$ -N ' + jobName + '\n\n' +
  
      "pushd " + jobDir + "\n\n") 

      file.write("setenv CCP4_SCR $TMPDIR\n\n")

      file.write('mrbump HKLIN ' + self.HKLIN + ' SEQIN ' + self.SEQIN +' HKLOUT ' + 'OUT.mtz  XYZOUT OUT.pdb << eof\n' +
      'LABIN ' + "F=" + self.LABIN["F"] + ' ' + "SIGF=" + self.LABIN["SIGF"] + ' ' + "FreeR_flag=" + self.LABIN["FreeR_flag"] + '\n' +
      'JOBID '+ str(jobID) + '_mrbump\n' +
      'MRPROGRAM '+mrbump_programs+'\n' +
      'LOCALFILE ' + ensemblePDB + ' CHAIN ALL RMS 1.2\n' +
    
      'SCOPSEARCH False\n' +
      'PQSSEARCH False\n' +
      'SSMSEARCH False\n' +
    
      'FAST False\n' +
      'DOFASTA False\n' +
    
      'MDLD False\n' +
      'MDLC False\n' +
      'MDLS False\n' +
      'MDLM False\n' +
      'MDLP False\n' +
    
      'FIXSG True\n' +
      'PJOBS 1\n' +
    
      'BUCC '  +       str(   self.BUCC )+'\n'+
      'BCYCLES '+     str(  self.BCYCLES )+'\n'+
      'ARPWARP '+       str(   self.ARPWARP )+'\n'+
      'ACYCLES '+     str(  self.ACYCLES )+'\n'+
      'SHELXE ' +      str(self.SHELXE)+'\n'+
      'SCYCLES '+       str( self.SCYCLES) +'\n' + 

      'CHECK False\n' +
      'LITE True\n' +
      'PICKLE False\n' +
      'UPDATE False\n' +
      'TRYALL True\n' +
      'USEACORN False\n' +
      'USEENSEM False\n' +
      'CLEAN False\n' +
      'DEBUG True\n') 

      if fixedPDB!="":
        file.write('FIXED_XYZIN ' + fixedPDB + ' IDEN ' + str(fixedIDEN) + '\n')

      file.write('END\n' +
      'eof\n\n' +
 
      'mkdir ' + os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx') +'\n' + 
      'mkdir ' + os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx') +'\n\n' + 

      'pushd ' + os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx') +'\n\n' +

      self.shelxClusterScript + ' ' + phaserPDB + ' ' + self.HKLIN  
                              + ' ' + str(matthews_coef.best_nmol) 
                              + ' ' + str(matthews_coef.best_solvent) 
                              + ' ' + str(self.MTZresolution) 
                              + " PHASER"
                              + " " + self.LABIN["F"] 
                              + " " + self.LABIN["SIGF"] 
                              + " " + self.LABIN["FreeR_flag"] + "\n\n" +
      'popd\n\n' +
      'cp ' + phaserPDB +' ' +os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx') +'\n' +
      'cp ' + phaserTFZPDB + ' '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx') +'\n\n' +
      'pushd ' + os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx') +'\n\n' +

      self.shelxClusterScript + ' ' + molrepPDB + ' ' + self.HKLIN  
                              + ' ' + str(matthews_coef.best_nmol) 
                              + ' ' + str(matthews_coef.best_solvent) 
                              + ' ' + str(self.MTZresolution) 
                              + " MOLREP"
                              + " " + self.LABIN["F"] 
                              + " " + self.LABIN["SIGF"] 
                              + " " + self.LABIN["FreeR_flag"] + "\n\n" +
      'popd\n\n' +
      'cp ' + molrepPDB+' ' + os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx') +'\n\n' +
    
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'data')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'logs')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'input')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'models')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'PDB_files')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'results')+'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'scratch') +'\n\n' +
      #'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'sequences')+'\n\n' +
#
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.ins')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.mtz') +'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.pro')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.ps')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.fcf')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.hkl')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.res')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.phs')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'phaser_shelx', 'orig.hat')+'\n\n' +
#
#
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.ins')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.mtz')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.pro')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.ps')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.fcf')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.hkl')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.res')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.phs')+'\n\n' + 
#      'rm -r '+os.path.join(jobDir, 'search_' + str(jobID) + '_mrbump', 'molrep_shelx', 'orig.hat')+'\n\n' +
#
#      'rm -r '+os.path.join(jobDir, 'OUT.pdb')+'\n\n' +
#      'rm -r '+os.path.join(jobDir, 'OUT.mtz')+'\n\n' +
      'popd\n' +  
      
      '\n\n')

      file.close()

      # Submit the job
      curDir=os.getcwd()
      os.chdir(jobDir)
      command_line='qsub -V %s' % sub_script
  
      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                    stdout = subprocess.PIPE)
  
      (child_stdout, child_stdin) = (p.stdout, p.stdin)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout.readline()
  
      qNumber=0
  
      while out:
         #sys.stdout.write(out)
         if self.QTYPE=="SGE":
            if "Your job" in out:
               qNumber=int(string.split(out)[2])
               self.qList.append(qNumber)
         out=child_stdout.readline()
  
      child_stdout.close()
  
      os.chdir(curDir)

   def makeResultsFile(self):
      """ Create final results summary file """
      
      resultsFile=os.path.join(self.RunDir, "Final_results.log")
  
      resfile=open(resultsFile, "w")

      for file in self.jobLogsList:
         if os.path.isfile(file):
            f=open(file, "r")
            line=f.readline()
            while line:
               if "SHELX>>>" in line:
                  resfile.write(line).replace("SHELX>>> ", "")
               line=f.readline()
      
            file.close()     
         else:
            resfile.write("Log file not found:\n   " + file + "\n")

      resfile.close()
      

if __name__ == "__main__":
   
   
   c=ClusterRun()
   c.QTYPE="SGE"
   c.qList=["2553845", "2553846", "2553847", "2553849"]
   c.monitorQueue(user="jac45")
   #c.monitorQueue()
   #list=c.getRunningJobList(user="jac45")
