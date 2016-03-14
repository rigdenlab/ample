#! /usr/bin/env python
#
#     Copyright (C) 2005 Ronan Keegan
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Application.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
#
# A script to run Matthews_coef 
# Ronan Keegan 30/07/07

import os, string, sys
import subprocess
import shlex
import re


class MattCoef:
   """ A class to run Matthews_coef. """
 
   def __init__(self):

      self.matt_EXE=os.path.join(os.environ["CCP4"], "bin", "matthews_coef")
      self.key=""

      self.best_prob=0.0
      self.best_nmol=0
      self.best_solvent=0.0

      self.debug=False

      self.CELL=""
      self.SYMM=""
      self.RESO=0.0

   def setCELL(self, cell_dimensions):
      self.CELL = cell_dimensions

   def setSYMM(self, symm):
      self.SYMM = symm

   def setRESO(self, resolution):
      self.RESO = resolution

   def runMC(self, target_MW, logfile, fixed_MW=0.0):
      """ Run Matthews_coef when mol weight includes a fixed component. 
          This looks at the probabilities while increasing the number of
          target components in the molecular weight. When a max is reached
          this value is taken as the number of target moleules in the a.s.u."""

      # Iterate, increasing the mw due to the target component until the best probability
      # is found
        
      # Calculate the new molecular weight
      mol_weight=target_MW + fixed_MW
   
      # Set the keywords
      self.key ="CELL %s\n" % self.CELL
      self.key+="SYMM %s\n" % self.SYMM
      self.key+="RESO %.3f\n" % self.RESO
      self.key+="MOLW %.3f\n" % mol_weight
      self.key+="AUTO\n"
      self.key+="END\n"
    
      command_line=self.matt_EXE

      # run Matthes_coef
      if os.name == "nt":
         process_args=shlex.split(command_line, posix=False)
         p = subprocess.Popen(process_args, shell="True", stdin = subprocess.PIPE,
                       stdout = subprocess.PIPE)
      else:
         process_args = shlex.split(command_line)
         p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                           stdout = subprocess.PIPE)

      (o, i) = (p.stdout, p.stdin)
      

      i.write(self.key)
      i.close()
    
      if os.path.isfile(logfile):
         log=open(logfile, "a")
      else:
         log=open(logfile, "w")
   
      line=o.readline()
      while line:
         if self.debug:
            sys.stdout.write(line)
   
         pattern = re.compile('^\s*(\d)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
         result = re.search(pattern, line)
         if result:
            split = string.split(line)
            nmol  = int(split[0])
            prob  = float(split[4])
            solvent  = float(split[2])
            if prob > self.best_prob:
               self.best_prob = prob
               self.best_nmol = nmol
               self.best_solvent = solvent
       
         line=o.readline()
      
      log.close()
      o.close()
    
if __name__ == "__main__":

   mc=MattCoef()

   # The 1l2w example (8 * protein "A" + 4 * protein "B")
   mc.setCELL("72.8400   73.3500   74.2600  103.4000  109.2000  107.4000")
   mc.setSYMM("P1")
   mc.setRESO(1.962)

   mc.runMC(13600, "matt_coef.log")

   print "Nmol:", mc.best_nmol, "Probability:", mc.best_prob, "Solvent content:", mc.best_solvent
