#!/usr/bin/env python

# Class for submitting ample modelling and MR/Refine/Build components to a cluster
# queuing system.
#
# Ronan Keegan 25/10/2011

import sys, os, string
import shlex, subprocess
import shutil


class Shelx:


   def __init__(self):
      self.pdbFile=""
      self.mtzFile=""
      self.phaserBestCC=0.0
      self.molrepBestCC=0.0
      self.phaserShelxPDB="Not set"
      self.molrepShelxPDB="Not set"

   def checkFileExist(self, filename, program, type):
      """ Check to see if the PDB file has been outout from Refmac in MrBUMP """

      if os.path.isfile(filename):
         status = 0
      else:
         status = 1
         sys.stdout.write("Warning: no file output from "+ program +"  for this job\n")
      
      return status

   def runMtz2hkl(self, rfree_flag=""):
      """ Run mtz2hkl to convert the MTZ file for Shelx """

      if rfree_flag != "":
         command_line='mtz2hkl -2 orig -r %s' % rfree_flag
      else:
         command_line='mtz2hkl -2 orig'

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      while out:
         sys.stdout.write(out)
         out=child_stdout_and_stderr.readline()
  
      child_stdout_and_stderr.close()

      status=self.checkFileExist("orig.hkl", "mtz2hkl", "hkl")
      if status==1:
         sys.stdout.write("Error: No output file from mtz2hkl\n")


   def runShelxpro(self, nmol):
      """ Run shelxpro to convert the PDB file for Shelx """

      command_line='shelxpro orig'

      keywords='I\n\norig.ins\norig.pdb\n\n\n' + str(nmol) + '\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nQ\n\n\nQ\n\n\n'

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.write(keywords)
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      while out:
         sys.stdout.write(out)
         out=child_stdout_and_stderr.readline()
  
      child_stdout_and_stderr.close()

      status=self.checkFileExist("orig.ins", "shelxpro", "ins")
      if status==1:
         sys.stdout.write("Error: No output file from shelxpro\n")


   def runShelxl(self):
      """ Run shelxl """

      command_line='shelxl orig'

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      while out:
         sys.stdout.write(out)
         out=child_stdout_and_stderr.readline()
  
      child_stdout_and_stderr.close()

      status=self.checkFileExist("orig.res", "shelxl", "res")
      if status==1:
         sys.stdout.write("Error: No output file from shelxl\n")



   def runShelxe(self, solvent, resolution, mrProgram):
      """ Run shelxe """

      if resolution<=2.0:
         free_lunch=" -e1.0 -l5"
      else:
         free_lunch=""

      frac_solvent=solvent/100.0

      command_line='shelxe  orig -a15 -q -s' + str(frac_solvent) + ' -h -z -b -l5' + free_lunch

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      while out:
         sys.stdout.write(out)

         if "CC for partial structure" in out:
            try:
               CCvalue=float(string.split(out)[-2])
            except:
               CCvalue=0.0
            if mrProgram.upper()=="PHASER": 
              if CCvalue > self.phaserBestCC:
                 self.phaserBestCC=CCvalue
                 self.phaserShelxPDB=os.path.join(os.getcwd(), "orig.pdb")
            if mrProgram.upper()=="MOLREP": 
              if CCvalue > self.molrepBestCC:
                 self.molrepBestCC=CCvalue
                 self.molrepShelxPDB=os.path.join(os.getcwd(), "orig.pdb")

         out=child_stdout_and_stderr.readline()
 
      child_stdout_and_stderr.close()

      status=self.checkFileExist("orig.phs", "shelxe", "phs")
      if status==1:
         sys.stdout.write("Error: No output file from shelxe\n")

   def results(self, mrProgram):
      """  Compile a results log for this Shelx job """

      sys.stdout.write("SHELX>>> \n")
      sys.stdout.write("SHELX>>> ########################################################################\n")
      sys.stdout.write("SHELX>>> Shelx results for PDB: " + self.pdbFile + "\n")
      sys.stdout.write("SHELX>>> \n")
      if mrProgram.upper() == "PHASER":
         sys.stdout.write("SHELX>>>    Best CC score for Phaser: " + str(self.phaserBestCC) + "\n")
         sys.stdout.write("SHELX>>>    Output PDB file: " + self.phaserShelxPDB + "\n")
      if mrProgram.upper() == "MOLREP":
         sys.stdout.write("SHELX>>>    Best CC score for Molrep: " + str(self.molrepBestCC) + "\n")
         sys.stdout.write("SHELX>>>    Output PDB file: " + self.molrepShelxPDB + "\n")
      sys.stdout.write("SHELX>>> \n")
      sys.stdout.write("SHELX>>> ########################################################################\n")
      sys.stdout.write("SHELX>>> \n")


if __name__ == "__main__":

   
   if  len(sys.argv) == 7:
      pdbFile=sys.argv[1]
      mtzFile=sys.argv[2]
      nmol=int(sys.argv[3])
      solvent=float(sys.argv[4])
      resolution=float(sys.argv[5])
      mrProgram=sys.argv[6]
      FreeR_flag=""
   elif  len(sys.argv) == 8:
      pdbFile=sys.argv[1]
      mtzFile=sys.argv[2]
      nmol=int(sys.argv[3])
      solvent=float(sys.argv[4])
      resolution=float(sys.argv[5])
      mrProgram=sys.argv[6]
      FreeR_flag=sys.argv[7]
   else:
      sys.stdout.write("Usage: python shelx_cluster.py <pdb> <mtz> <nmol> <solvent> <mtz resolution> <MR program> <FreeR_flag (optional)>\n")
      sys.stdout.write("\n")
      sys.stdout.write("       e.g. python shelx_cluster.py foo.pdb foo.mtz 2 45.33 1.543 PHASER FreeR_flag\n")
      sys.stdout.write("\n")
      sys.exit()
     
   SHELX=Shelx()

   status=SHELX.checkFileExist(pdbFile, "mrbump", "pdb")
   if status==0:
      shutil.copy(pdbFile, "orig.pdb")
      SHELX.pdbFile="orig.pdb"
   else:
      sys.exit(1)

   status=SHELX.checkFileExist(mtzFile, "mtzfile", "mtz")
   if status==0:
      shutil.copy(mtzFile, "orig.mtz")
      SHELX.mtzFile="orig.mtz"
   else:
      sys.exit(1)

   SHELX.runMtz2hkl(rfree_flag=FreeR_flag)
   SHELX.runShelxpro(1)
   SHELX.runShelxl()
   SHELX.runShelxe(solvent, resolution, mrProgram)
   SHELX.results(mrProgram)

