#!/usr/bin/env python
#
#     Copyright (C) 2012 Ronan Keegan 
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Application.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
# Class for running Shelxe in Ample for c-alpha tracing
#
# Ronan Keegan 04/05/2012

import sys, os, string
import shlex, subprocess
import shutil

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# Test for environment variables and required executables

if not "CCP4" in sorted(os.environ.keys()):
    raise RuntimeError('CCP4 not found')
if not which("mtz2various"):
    raise RuntimeError('mtz2various not found')
if not which("shelxe"):
    raise RuntimeError('shelxe not found')

class Shelx:

   def __init__(self):
      self.pdbinFile=""
      self.mtzinFile=""
      self.hklinFile=""
      self.pdboutFile=""
      self.phsoutFile=""

      self.phaserBestCC=0.0
      self.molrepBestCC=0.0
      self.phaserShelxPDB="Not set"
      self.molrepShelxPDB="Not set"
      self.shelxScriptFile="shelx-script.sh"
      self.script=""
      self.PDBspacegroup=""
      self.PDBCRYST1Line=""

      self.labin=dict([])

      try:
         self.debug=eval(os.environ['AMPLE_DEBUG'])
      except:
         self.debug=False

      self.mtz2variousLogfile="mtz2various.log"
      self.shelxeLogfile="shelxe.log"
     
   def setDebug(self, debugValue):
      self.debug=debugValue

   def setShelxScriptFile(self, filename):
      self.shelxScriptFile=filename

   def checkFileExist(self, filename, program, type):
      """ Check to see if the PDB file has been outout from Refmac in MrBUMP """

      if os.path.isfile(filename):
         status = 0
      else:
         status = 1
         sys.stdout.write("Warning: no file output from "+ program +"  for this job\n")
      
      return status

   def parsePDBfileSG(self, pdbfile):
      """ Check to see if the space group is P -1 """
   
      cryst1Line=""
      if os.path.isfile(pdbfile):
         f=open(pdbfile, "r")
         line=f.readline()
         while line:
            if "CRYST1" in line[0:6]:
               self.PDBspacegroup=string.join(string.split(line)[7:]).replace(" ","").upper()
               self.PDBCRYST1Line=line
               if self.PDBspacegroup == "P-1":
                  sys.stdout.write("Warning: Shelxe - changing input spacegroup in pda file from P -1 to P 1\n")
                  cryst1Line=line.replace("-1", "1")
            line=f.readline()
   
         f.close()
         if self.PDBspacegroup == "P-1":
            f=open(pdbfile, "r")
            lines=f.readlines()
            f.close()
            import random
            tempfile=os.path.join(os.environ["CCP4_SCR"], "shelxe-pdb_sgp1_" + str(random.random()).replace(".","") + ".tmp")
            o=open(tempfile,"w")
            for line in lines:
               if "CRYST1" not in line[0:6]:
                  o.write(line)
               else:
                  o.write(cryst1Line)
            o.close()
            os.remove(pdbfile)
            shutil.move(tempfile, pdbfile)
      else:
         sys.stdout.write("Warning: Shelxe pdb file correction - pdb file could not be found\n")
   
   def runMtz2Various(self, mtzinFile, fp="FP", sigfp="SIGFP", free="FreeR_flag"):
      """ Run mtz2various to convert the MTZ file for Shelxe """

      self.mtzinFile=mtzinFile
      self.hklinFile=os.path.splitext(mtzinFile)[0] + ".hkl"

      command_line='mtz2various HKLIN %s HKLOUT %s' % (self.mtzinFile, self.hklinFile)

      self.script += command_line + " << eof\n"

      key  = "LABIN FP=%s SIGFP=%s FREE=%s\n" % (fp, sigfp, free)
      key += "OUTPUT SHELX\n"
      key += "FSQUARED\n"
      key += "END\n"
   
      self.script += key + "eof\n\n"

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.write(key)
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      if self.debug: 
         log=open(self.mtz2variousLogfile, "w")

      while out:
         sys.stdout.write(out)
         out=child_stdout_and_stderr.readline()
         if self.debug:
            log.write(out)
  
      child_stdout_and_stderr.close()
      if self.debug:
         log.close()

      status=self.checkFileExist(self.hklinFile, "mtz2various", "hkl")
      if status==1:
         sys.stdout.write("Error: No output file from mtz2various\n")


   def runMtz2hkl(self, rfree_flag=""):
      """ Run mtz2hkl to convert the MTZ file for Shelx """

      if rfree_flag != "":
         command_line='mtz2hkl -2 orig -r %s' % rfree_flag
      else:
         command_line='mtz2hkl -2 orig'

      self.script += command_line + "\n"

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

   def runShelxe(self, solvent, resolution, mrProgram, pdbinFile, pdboutFile="shelxe-output.pdb", phsoutFile="shelxe-output.phs", traceCycles=20):
      """ Run shelxe """

      self.pdbinFile=pdbinFile
      self.pdboutFile=pdboutFile
      self.phsoutFile=phsoutFile

      if resolution<=2.0:
         free_lunch=" -e1.0 -l5"
      else:
         free_lunch=""

      frac_solvent=solvent/100.0

      # Set up input files for shelxe
      if os.path.isfile("shelxe-input.hkl"):
         os.remove("shelxe-input.hkl")
      shutil.copyfile(self.hklinFile, "shelxe-input.hkl")
      self.script += "cp %s shelxe-input.hkl\n" % self.hklinFile
      if os.path.isfile("shelxe-input.pda"):
         os.remove("shelxe-input.pda")
      shutil.copyfile(self.pdbinFile, "shelxe-input.pda")
      #Check the input PDA file for SG P -1
      self.parsePDBfileSG("shelxe-input.pda")
      self.script += "cp %s shelxe-input.pda\n\n" % self.pdbinFile

      command_line='shelxe shelxe-input.pda -a' + str(traceCycles) + ' -q -s' + str(frac_solvent) + ' -h -z -b -l5' + free_lunch
      self.script += command_line + "\n\n"

      file=open(self.shelxScriptFile, "w")
      file.write(self.script) 
      file.close()

      process_args = shlex.split(command_line)
      p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
  
      (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
  
      # Write the keyword input
      child_stdin.close()
  
      # Watch the output for successful termination
      out=child_stdout_and_stderr.readline()
  
      if self.debug:
         log=open(self.shelxeLogfile, "w")
      while out:
         sys.stdout.write(out)
         if self.debug:
            log.write(out)
         out=child_stdout_and_stderr.readline()
 
      child_stdout_and_stderr.close()
      if self.debug:
         log.close()

      status=self.checkFileExist("shelxe-input.pdb", "shelxe", "pdb")
      if status==1:
         sys.stdout.write("Error: No output pdb file from shelxe\n")
      else:
         shutil.move("shelxe-input.pdb", self.pdboutFile)
         script  = "mv shelxe-input.pdb %s\n" % self.pdboutFile

      status=self.checkFileExist("shelxe-input.phs", "shelxe", "phs")
      if status==1:
         sys.stdout.write("Error: No output phs file from shelxe\n")
      else:
         shutil.move("shelxe-input.phs", self.phsoutFile)
         script += "mv shelxe-input.phs %s\n\n" % self.phsoutFile

      scriptFile=open(self.shelxScriptFile, "a")
      scriptFile.write(script)
      scriptFile.close()
 
      # Set the best scores
      if os.path.isfile(self.pdboutFile):
         score=0.0
         f=open(self.pdboutFile, "r")
         line=f.readline()
         try:
            score=float(string.split(line)[6].replace("%",""))
         except: 
            sys.stdout.write("Warning (shelxe_trace): Can't find Shelxe output score in output pdb\n")
         f.close()
         if mrProgram.upper() == "PHASER": 
            self.phaserBestCC=score 
         if mrProgram.upper() == "MOLREP": 
            self.molrepBestCC=score 

         # Set SG back to P -1 if that what is what it was originally
         if self.PDBspacegroup=="P-1":
            f=open(self.pdboutFile, "r")
            lines=f.readlines()
            f.close()
            import random
            tempfile=os.path.join(os.environ["CCP4_SCR"], "shelxe-pdb_sgp-1_" + str(random.random()).replace(".","") + ".tmp")
            o=open(tempfile, "w")
            for line in lines:
               if "CRYST1" not in line:
                  o.write(line)
               else:
                  o.write(self.PDBCRYST1Line)
            o.close()
            os.remove(self.pdboutFile)
            shutil.move(tempfile, self.pdboutFile)


      if self.debug == False:
        for file in self.hklinFile, "shelxe-input.hkl", "shelxe-input.pdb", \
                    "shelxe-input.pha", "shelxe-input.hat", \
                    phsoutFile, "shelxe-input.pda": 
           if os.path.isfile(file):
              os.remove(file) 

   def results(self, mrProgram):
      """  Compile a results log for this Shelx job """

      sys.stdout.write("SHELX>>> \n")
      sys.stdout.write("SHELX>>> ########################################################################\n")
      sys.stdout.write("SHELX>>> Shelx results for PDB: " + self.pdbinFile + "\n")
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
      pdbinFile=sys.argv[1]
      mtzinFile=sys.argv[2]
      nmol=int(sys.argv[3])
      solvent=float(sys.argv[4])
      resolution=float(sys.argv[5])
      mrProgram=sys.argv[6]
      FP="FP"
      SIGFP="SIGFP"
      FREE="FREE"
   elif  len(sys.argv) == 10:
      pdbinFile=sys.argv[1]
      mtzinFile=sys.argv[2]
      nmol=int(sys.argv[3])
      solvent=float(sys.argv[4])
      resolution=float(sys.argv[5])
      mrProgram=sys.argv[6]
      FP=sys.argv[7]
      SIGFP=sys.argv[8]
      FREE=sys.argv[9]
   else:
      sys.stdout.write("Usage: python shelxe_trace.py <pdb> <mtz> <nmol> <solvent> <mtz resolution> <MR program> [<FP> <SIGFP> <FREE>]\n")
      sys.stdout.write("\n")
      sys.stdout.write("       e.g. python shelxe_trace.py foo.pdb foo.mtz 2 45.33 1.543 PHASER FP SIGFP FREE\n")
      sys.stdout.write("\n")
      sys.exit()
     
   SHELX=Shelx()

   status=SHELX.checkFileExist(pdbinFile, "mrbump", "pdb")
   if status==1:
      sys.exit(1)

   status=SHELX.checkFileExist(mtzinFile, "mtzfile", "mtz")
   if status==1:
      sys.exit(1)

   SHELX.runMtz2Various(mtzinFile, fp=FP, sigfp=SIGFP, free=FREE)
   SHELX.runShelxe(solvent, resolution, mrProgram, pdbinFile, traceCycles=15)
   SHELX.results(mrProgram)

