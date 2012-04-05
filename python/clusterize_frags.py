#!/usr/bin/env python

# for submitting ample modelling and MR/Refine/Build components to a cluster
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



###############
def make_frag_script(cmd, jobDir, jobID, jobName,  scriptname):        
      os.chdir(jobDir)
      os.mkdir( os.path.join(jobDir, "logs") )

      script=open(scriptname, "w")
      script.write('#!/bin/sh\n'
         '#$ -j y\n' +
         '#$ -cwd\n' +
         '#$ -w e\n' +
         '#$ -V\n' +
         '#$ -S /bin/sh\n'+
         '#$ -o ' + os.path.join(jobDir, "logs", "job_" + str(jobID) + '.log') + '\n' +
         '#$ -N A' + jobName + '\n\n' +
          "pushd " + jobDir + "\n\n") 

      script.write(cmd) 

      curDir=os.getcwd()
      os.chdir(jobDir)


      ##submit
      command_line='qsub '+ scriptname
  
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
            if "NNMAKE DONE" in out:
               qNumber=int(string.split(out)[2])
               self.qList.append(qNumber)
         out=child_stdout.readline()
  
      child_stdout.close()  
      os.chdir(curDir)

if __name__=="__main__":
   

 cmd =  '/home/jac45/rosetta/rosetta_fragments/nnmake/make_fragments.pl 1EJG_.fasta -noprof -nojufo -nosam -v -id 1EJG_'
 jobdir = '/home/jac45/1EJG'
 jobid= 'test'
 scriptname = '/home/jac45/1EJG/1EJG_run.sub'
 jobname='jobname'
 make_frag_script(cmd, jobdir, jobid, jobname, scriptname)
  
