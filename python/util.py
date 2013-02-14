'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''
import os
import sys
import re

import unittest

import MTZParse

# get a program test for existsnce
def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

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
#########
def check_for_exe(exename, varname):
    exepath = ''
    print 'looking for', exename
    if not varname:
        print 'no ' + exename + ' given on the command line, looking in the PATH'
        print which(exename)
        if not which(exename):
            print 'You need to give the path for ' + exename
            sys.exit()
        else:

            exepath = which(exename)

    else:
        exepath = varname
        print 'using here', exepath




    if not os.path.exists(exepath):
        print 'You need to give the path for ' + exename + ', executable in the PATH dosnt exist'
        sys.exit()
    else:
        return exepath
########

def get_mtz_flags( mtzfile ):
    """
    Return the SIGF, F and FREER column labels in an MTZ file
    """
    
    # Instantiate the class
    md=MTZParse.Mtzdump()
    
    # Set the MTZ file and the output log file
    md.setHKLIN( mtzfile )
    md.setMTZdumpLogfile("./")
    
    # Run mtzdump
    md.go()
    
    # Extract column data
    md.getColumnData()
    
    F =  md.colLabels['F']
    SIGF =  md.colLabels['SIGF']
    FREE =  md.colLabels['FREE']

    return (F, SIGF, FREE)
    
    
##End get_mtz_flags

def XXXget_mtz_flags(mtz):
    """
    Old version of get_flags
    """
    sigf = 'SIGF='
    FP = 'F='
    free = 'FreeR_flag=Unassigned'
    
    path = os.getcwd()
    os.system('mtzdmp ' + mtz + ' >mtzdmp_out')
    mtz_out = open(path + '/mtzdmp_out')
    flags = False
    while flags == False:
        line = mtz_out.readline()
        get_stat1 = re.compile('Column Labels')
        result_stat1 = get_stat1.search(line)
        if result_stat1:
           #print line
           next = mtz_out.readline()
           next = mtz_out.readline()
        if  re.search('Column Types', line):
           lab = mtz_out.readline()
           lab = mtz_out.readline() 
    
           flags = True
    #  print next
    #  print lab
    
    
    
    
    lab_list=re.split('\s*', lab)
    name_list=re.split('\s*', next)
    
    sigf = lab_list.index('Q')
    sigf = name_list[sigf]
    
    FP = lab_list.index('F')
    FP = name_list[FP]
    
    free = lab_list.index('I')
    free = name_list[free]
    
    # print sigf, FP, free
    
    #jmht return next, sigf, FP, free
    return next, FP, sigf, free


def get_psipred_prediction(psipred):
    string = ''
    for line in open(psipred):
        get_stat1 = re.compile('^Pred\:\s*(\w*)')
        result_stat1 = get_stat1.match(line)
        if result_stat1:
            stat1_get = re.split(get_stat1, line)
            # print stat1_get[1]
            string = string + stat1_get[1]


    C = 0
    H = 0
    E = 0
    length = len(string)

    for c in string:
        if c == 'C':
            C = C + 1
        if c == 'H':
            H = H + 1
        if c == 'E':
            E = E + 1


    H_percent = float(H) / length * 100
    E_percent = float(E) / length * 100

    if H > 0 and E > 0:
        print  'Your protein is predicted to be mixed alpha beta, your chances of success are intermediate'
    if H == 0 and E > 0:
        print  'Your protein is predicted to be all beta, your chances of success are low'
    if H > 0 and E == 0:
        print  'Your protein is predicted to be all alpha, your chances of success are high'
    if  H == 0 and E == 0:
        print  'Your protein is has no predicted secondary structure, your chances of success are low'
        


class TestCell(unittest.TestCase):
    """
    Unit test
    """
    
    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """
                
        thisdir = os.getcwd()
        self.ampledir = os.path.abspath( thisdir+os.sep+"..")
        self.testdir = self.ampledir + os.sep + "tests"
        self.testfilesdir = self.testdir + os.sep + "testfiles"
    
    def testGetMtzFlags(self):
        """Get MTZ flags"""
        
        
        mtz = self.ampledir + os.sep + "examples" + os.sep + "toxd-example" + os.sep + "1dtx.mtz"
        
        t_flag_F, t_flag_SIGF, t_flag_FREE = get_mtz_flags( mtz )
        
        f = 'FP'
        sigf = 'SIGFP'
        free = 'FreeR_flag'
        
        self.assertEqual(f, t_flag_F, "Correct F")
        self.assertEqual(sigf, t_flag_SIGF, "Correct SIGF")
        self.assertEqual(free, t_flag_FREE, "Correct FREE")
        
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
