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

import logging
import os
import shlex
import subprocess
import sys
import time


class Mtzdump:

   """This is a wrapper for the program mtzdump which will print the
   header of an MTZ file. Methods are provided for extracting certain
   quantities from the header"""
 
   def __init__(self):
     self.mtzdumpEXE=os.path.join(os.environ["CBIN"], "mtzdump")
     self.programParameters = 'END\n'
     self.logfile=""
     self.hklin=""
     self.colLabels=dict([])
     try:
        self.debug=eval(os.environ['MRBUMP_DEBUG'])
     except:
        self.debug=False

     self.logger = logging.getLogger()
   
   def setMTZdumpLogfile(self, tdir):
     """ Set the name of the mtzdump logfile. Try to make it unique to the user
         by appending their username to the start of the file name. """
     if os.name == "nt":
        self.logfile = os.path.join(tdir, os.environ["USERNAME"] + '_junk_mtzdump.log')
     else:
        self.logfile = os.path.join(tdir, os.environ["USER"] + '_junk_mtzdump.log')
    
   def setProgramParameters(self,inputline):
     self.programParameters = inputline + '\n' + self.programParameters
 
   def setHKLIN(self, hklin):
     self.hklin = hklin
 
   def go(self):
      """ Run Mtzdump """
 
      # Set the command line
      command_line = self.mtzdumpEXE + ' HKLIN ' + self.hklin
         
#      if self.debug:
#         sys.stdout.write("\n")
#         sys.stdout.write("======================\n")
#         sys.stdout.write("MTZDUMP command line:\n")
#         sys.stdout.write("======================\n")
#         sys.stdout.write(command_line + "\n")
#         sys.stdout.write("\n")

      self.logger.debug("Executing MTZDUMP with command-line: {0}".format( command_line ) )
 
      # Launch
      if os.name == "nt":
         process_args = shlex.split(command_line, posix=False)
         p = subprocess.Popen(process_args, shell="True", stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE)
      else:
         process_args = shlex.split(command_line)
         p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE)
 
      (child_stdout, child_stdin) = (p.stdout, p.stdin)
 
      child_stdin.write(self.programParameters)
      child_stdin.close()         
 
      # Capture any error from the alignment job and print it to screen if debug mode is on
      log=open(self.logfile, "w")
      line=child_stdout.readline()
      while line:
         log.write(line)
         if self.debug:
            sys.stdout.write(line)
         line=child_stdout.readline()
      child_stdout.close()
      log.close()
 
      count=0
      # Wait for the logfile to be written
      while os.path.isfile(self.logfile) == False:
         time.sleep(1) 
         count=count+1
         if count > 30:
           sys.stdout.write("MTZdump Error: taking too long to write log file %s" % self.logfile)
           break
 
   def getCell(self):
     f = open(self.logfile,'r')
     logline = f.readline()
     # with "while logline" this loop terminates early
     # don't know why!
     while 1:
       if 'Cell Dimensions' in logline:
         logline = f.readline()
         logline = f.readline()
         cell = logline.split()
         f.close()
         return cell
       logline = f.readline()
     f.close()
 
   def getSymmetryNumber(self):
     """ Get the symmetry number. """
 
     f = open(self.logfile,'r')
     logline = f.readline()
     # with "while logline" this loop terminates early
     # don't know why!
     while 1:
       if 'Space group' in logline:
         l = logline.split()
         symmetryNumber = int(l[-1].rstrip(')'))
         f.close()
         return symmetryNumber
       logline = f.readline()
     f.close()
 
   def getSpacegroup(self):
     """ Get the space group """
 
     f = open(self.logfile,'r')
     logline = f.readline()
 
     while logline:
       if 'Space group' in logline:
         l = logline.split("'")
         spacegroup = l[-2].replace(" ","")
         f.close()
         return spacegroup
       logline = f.readline()
     f.close()
 
   def getResolution(self):
     """ Get the resolution of the target MTZ file. """
 
     f = open(self.logfile,'r')
     logline = f.readline()
     # with "while logline" this loop terminates early
     # don't know why!
     while 1:
       if 'Resolution Range' in logline:
         logline = f.readline()
         logline = f.readline()
         resline = logline.split()
         resolution = float(resline[-3].rstrip(')'))
         f.close()
         return resolution
       logline = f.readline()
     f.close()
 
   def getColumnData(self):
     """ A function to get the column labels and types from an MTZ file. Takes in a dictionary
     to store the labels and the their corresponding column types"""
 
     # Open the log file for reading
 
     f = open(self.logfile,'r')
     logline = f.readline()
 
     # Loop over the lines in the logfile and save the column data information
     col_labels = None 
     while logline:
       if 'Column Labels' in logline:
         logline = f.readline()
         logline = f.readline()
         col_labels = logline.strip().split()
       #if 'Column Types' in logline:
       #  logline = f.readline()
       #  logline = f.readline()
       #  ctypes = string.split(string.strip(logline))
         break
       logline = f.readline()

     if not col_labels:
         msg = """Error parsing MTZDUMP logfile, could not find "Column Labels" header.
Please check the logfile: {0} for any errors.""".format( self.logfile )
         self.logger.critical(msg)
         raise RuntimeError,msg

     for col in col_labels:
        if col == "F":
           self.colLabels["F"]="F"
        if col == "FP":
           self.colLabels["F"]="FP"
        if col == "SIGF":
           self.colLabels["SIGF"]="SIGF"
        if col == "SIGFP":
           self.colLabels["SIGF"]="SIGFP"
        if col == "FREE":
           self.colLabels["FREE"]="FREE"
        if col == "FreeR_flag":
           self.colLabels["FREE"]="FreeR_flag"
             
   
     # Fill the dictionary with the column details
 
     #for i in range(len(clabels)):
        #col_data[clabels[i]]=ctypes[i]
 
     f.close()
     return col_labels
 
 
if __name__ == "__main__":
   """ A test example or can be used to run as standalone. """

   # Check the input arguments
   if len(sys.argv) != 3:
      print "Usage: python MTZParse.py <HKLIN> <info type>"
      print "       HKLIN     - input MTZ file"
      print "       info type - information requested"
      print "          -> cell"
      print "          -> spacegroup"
      print "          -> resolution"
      print "          -> col_data"
      sys.exit()

   # Set information type
   itype=(sys.argv[2]).lower()

   if "cell" not in itype and "spacegroup" not in itype and "resolution" not in itype and "col_data" not in itype:
      print "'%s' is not a valid argument for <info type>" % sys.argv[2]
      print "Valid arguments are:" 
      print "          -> cell"
      print "          -> spacegroup"
      print "          -> resolution"
      print "          -> col_data"
      sys.exit()

   # Instantiate the class
   md=Mtzdump()
   
   # Set the MTZ file and the output log file
   md.setHKLIN(sys.argv[1])
   md.setMTZdumpLogfile("./")

   # Run mtzdump
   md.go()
  
   if itype == "cell": 
      print md.getCell()
   
   if itype == "spacegroup": 
      print md.getSpacegroup()
 
   if itype == "resolution": 
      print md.getResolution()
 
   if itype == "col_data": 
      print md.getColumnData()
 

