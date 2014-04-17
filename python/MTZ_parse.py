#! /usr/bin/env ccp4-python
#
#     Copyright (C) 2005 Ronan Keegan
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Application.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
#

import os
import sys
import string
import subprocess
import shlex

class MTZ_parse:
    """ A class to parse and store the output from mtzdmp. """

    def __init__(self):
        self.mtzdmp_EXE  = os.path.join(os.environ["CBIN"], "mtzdmp")
        self.mtzdump_EXE = os.path.join(os.environ["CBIN"], "mtzdump")
        self.MTZ_file = ''
        self.no_of_columns = 0
        self.no_of_reflections = 0

        self.cell_dimensions=dict([])
        self.resolution = 0.0
        self.space_group = ''
        self.log=[]
        self.col_labels=[]
        self.col_types=[]
        self.F=""
        self.SIGF=""
        self.FreeR_flag=""
        self.FreeR_valid=False
        try:
            self.debug=eval(os.environ['MRBUMP_DEBUG'])
        except:
            self.debug=False

    def checkRFREE(self, FreeR_flag=None):
        """Check the RFREE flag is valid"""
        
        # Check the given flag or use the one found by parsing the file
        if FreeR_flag is None:
            FreeR_flag = self.FreeR_flag
        
        # Need something to check...
        if FreeR_flag is None or len(FreeR_flag) == 0:
            return False
        
        # Make sure its there
        if FreeR_flag not in self.col_labels:
            return False
        
        # We search for the file statistics and check that the min and max rfree values
        # are different
        
        # Find the summary
        stats=-1
        assert len(self.log)
        for i, line in enumerate(self.log):
            if 'OVERALL FILE STATISTICS' in line:
                stats = i
                break
        
        if stats == -1:
            sys.stdout.write("Warning: could not find OVERALL FILE STATISTICS section in input MTZ file\n")
            sys.stdout.write("\n")
            return False
        #
        # Assume the data we want starts 7 lines in
        start = stats + 7
        try:
            for line in self.log[start:]:
                tokens = line.split()
                if tokens[11] == FreeR_flag:
                    minf = float(tokens[2])
                    maxf = float(tokens[3])
                    if minf != maxf:
                        return True
                    break
        except Exception,e:
            sys.stdout.write("Warning: error checking FreeR_valid in input MTZ file\n")
            sys.stdout.write("\n")
            
        return False

    def setDEBUG(self, flag):
        self.debug=flag

    def setMTZfile(self, MTZ_file):
        self.MTZ_file =  MTZ_file

    def setNoofCols(self, no_of_columns):
        self.no_of_columns = no_of_columns

    def setNoofReflections(self, no_of_reflections):
        self.no_of_reflections = no_of_reflections

    def setCellDimensions(self, cell_dimensions_array):
        self.cell_dimensions['a'] = float(cell_dimensions_array[0])
        self.cell_dimensions['b'] = float(cell_dimensions_array[1])
        self.cell_dimensions['c'] = float(cell_dimensions_array[2])
        self.cell_dimensions['alpha'] = float(cell_dimensions_array[3])
        self.cell_dimensions['beta'] = float(cell_dimensions_array[4])
        self.cell_dimensions['gamma'] = float(cell_dimensions_array[5])

    def setResolution(self, resolution_array):
        self.resolution = float(resolution_array[-3])

    def setSpaceGroup(self, space_group_array):
        self.space_group = space_group_array[1].replace(' ','')

    def getColumnData(self):
        """ A function to get the column labels and types from an MTZ file. Takes in a dictionary
        to store the labels and the their corresponding column types"""

        # Assumption is FREER data is duff
        self.FreeR_valid = False

        # Loop over the lines in the logfile and save the column data information
        count=0
        for line in self.log:
            if 'Column Labels' in line:
                self.col_labels = string.split(string.strip(self.log[count+2]))
                self.col_types = string.split(string.strip(self.log[count+6]))
                break
            count=count+1

        # Search for F and SIGF
        count=0
        for i in self.col_types:
            if i == "F":
                self.F=self.col_labels[count]
                self.SIGF="SIG"+self.F
                if self.SIGF not in self.col_labels:
                    sys.stdout.write("Error: SIGF label %s not found for %s in input MTZ file\n" % (self.SIGF, self.F))
                    sys.stdout.write("\n")
                    sys.exit()
                break
            count=count+1

        # Search for FreeR_flag
        count=0
        if self.col_types.count("I") > 1:
            sys.stdout.write("Warning: More than 1 free set found in input MTZ file\n")
            sys.stdout.write("\n")
        for i in self.col_types:
            if i == "I":
                FreeR_flag=self.col_labels[count]
                # Check this Rfree is valid
                if self.checkRFREE(FreeR_flag=FreeR_flag):
                    self.FreeR_flag=FreeR_flag
                    break
            count=count+1
        
        return

    def run_mtzdmp(self, MTZ_file_path):
        """ Run MTZDump on the input MTZ file to capture various information """

        if self.debug:
            sys.stdout.write("\n")
            sys.stdout.write("Preparation: Running MTZDump on the input MTZ file")
            sys.stdout.write("\n")

        self.setMTZfile(os.path.split(MTZ_file_path)[-1])
        
        # Need to clear log each time we run or we are just processing the last file
        self.log = []

        # Set the command line
        command_line = self.mtzdump_EXE + " HKLIN " + MTZ_file_path
        
        if self.debug:
            sys.stdout.write("\n")
            sys.stdout.write("======================\n")
            sys.stdout.write("MTZDUMP command line:\n")
            sys.stdout.write("======================\n")
            sys.stdout.write(command_line + "\n")
            sys.stdout.write("\n")

        # Launch Mafft
        if os.name == "nt":
            process_args = shlex.split(command_line, posix=False)
            p = subprocess.Popen(process_args, shell="True", stdin = subprocess.PIPE,
                                   stdout = subprocess.PIPE)
        else:
            process_args = shlex.split(command_line)
            p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                   stdout = subprocess.PIPE)

        (child_stdout, child_stdin) = (p.stdout, p.stdin)

        child_stdin.write("END\n")
        child_stdin.close()

        # Capture various pieces of useful information
        line = child_stdout.readline()
        while line:
            if "* Number of Columns" in line:
                self.setNoofCols(int(string.split(line)[-1]))
            if "* Number of Reflections" in line:
                self.setNoofReflections(int(string.split(line)[-1]))
            if "* Cell Dimensions :" in line:
                line = child_stdout.readline()
                line = child_stdout.readline()
                self.setCellDimensions(string.split(line))
            if "Resolution Range :" in line:
                line = child_stdout.readline()
                line = child_stdout.readline()
                self.setResolution(string.split(line))
            if "* Space group =" in line:
                self.setSpaceGroup(string.split(line,"'"))
            if self.debug:
                sys.stdout.write(line)

            self.log.append(string.strip(line))
            line = child_stdout.readline()

        child_stdout.close()
        self.getColumnData()
        
        return

if __name__ == '__main__':

    if len(sys.argv)!=2:
        sys.stdout.write("Usage: python MTZ_parse.py <mtzfile>\n")
        sys.exit()
    mtzfile=sys.argv[1]

    mtz = MTZ_parse()
    mtz.run_mtzdmp(mtzfile)
    print "cell_dimensions ",mtz.cell_dimensions
    print mtz.resolution
    print mtz.space_group
    print "col_labels ",mtz.col_labels
    print "col_labels ",mtz.col_types
    print "F ",mtz.F
    print "SIGF ",mtz.SIGF
    print "FreeR_flag ",mtz.FreeR_flag
    print "FreeR_valid ",mtz.FreeR_valid
