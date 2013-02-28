'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# Python modules
import logging
import os
import re
import sys
import unittest

# Our modules
import MTZParse

references = """AMPLE: To be added
SHELX: is used: "A short history of SHELX". Sheldrick, G.M. (2008). Acta Cryst. A64, 112-122

SCWRL: G. G. Krivov, M. V. Shapovalov, and R. L. Dunbrack, Jr. Improved prediction of protein
side-chain conformations with SCWRL4. Proteins (2009).

Theseus: THESEUS: Maximum likelihood superpositioning and analysis of macromolecular structures.
Theobald, Douglas L. & Wuttke, Deborah S. (2006b) Bioinformatics 22(17):2171-2172 [Open Access]
Supplementary Materials for Theobald and Wuttke 2006b.

MrBUMP: R.M.Keegan and M.D.Winn (2007) Acta Cryst. D63, 447-457

CCP4: Collaborative Computational Project, Number 4. (1994), The CCP4 Suite: Programs
for Protein Crystallography. Acta Cryst. D50, 760-763\n

MOLREP: A.A.Vagin & A.Teplyakov (1997) J. Appl. Cryst. 30, 1022-1025\n

PHASER: McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D.,
Storoni, L.C. & Read, R.J. (2007)
Phaser crystallographic software J. Appl. Cryst. 40, 658-674

REFMAC: G.N. Murshudov, A.A.Vagin and E.J.Dodson, (1997) Refinement of Macromolecular
Structures by the Maximum-Likelihood Method. Acta Cryst. D53, 240-255

SPICKER: Y. Zhang, J. Skolnick, SPICKER: Approach to clustering protein structures for
near-native model selection, 
Journal of Computational Chemistry, 2004 25: 865-871

MaxCluster: http://www.sbg.bio.ic.ac.uk/maxcluster/"""

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
    logger = logging.getLogger()
    exepath = ''
    logger.debug('looking for: {}'.format(exename) )
    if not varname:
        logger.debug( 'no {} given on the command line, looking in the PATH'.format(exename) )
        logger.debug( "{}".format(which(exename)) )
        if not which(exename):
            logger.critical('You need to give the path for: {}'.format(exename))
            sys.exit()
        else:

            exepath = which(exename)
    else:
        exepath = varname
        logger.debug( 'using here {}'.format(exepath) )

    if not os.path.exists(exepath):
        logger.critical( 'You need to give the path for {}, executable in the PATH dosnt exist'.format(exename) )
        sys.exit()
    else:
        return exepath
########

def make_workdir(work_dir, rootname='ROSETTA_MR_'):
    """
    Make a work directory rooted at work_dir and return its path
    """
    print 'Making a Run Directory  -checking for previous runs'
    run_inc = 0
    run_making_done = False
    while not run_making_done:
        if not os.path.exists(work_dir + os.sep + rootname + str(run_inc)):
            run_making_done = True
            os.mkdir(work_dir + os.sep +rootname + str(run_inc))
        run_inc += 1
    work_dir = work_dir + os.sep + rootname + str(run_inc - 1)
    return work_dir

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
        

def setup_logging():
    """
    Set up the various log files/console logging
    and return the logger
    """
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # create file handler and set level to debug
    fl = logging.FileHandler("JENS_AMPLE.log")
    fl.setLevel(logging.DEBUG)
    
    # create formatter for fl
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # add formatter to fl
    fl.setFormatter(formatter)
    
    # add fl to logger
    logger.addHandler(fl)
    
    # Now create console logger for outputting stuff
    # create file handler and set level to debug
    cl = logging.StreamHandler(stream=sys.stdout)
    cl.setLevel(logging.INFO)
    
    # create formatter for fl
    formatter = logging.Formatter('%(message)s')
    
    # add formatter to fl
    cl.setFormatter(formatter)
    
    # add fl to logger
    logger.addHandler(cl)
    
    return logger
        


class TestUtil(unittest.TestCase):
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
        
    def testGetMtzFlagsNoFree(self):
        """Get MTZ flags when there are no free flags"""
        
        mtz = self.testfilesdir + os.sep + "2uui_sigmaa.mtz"
        self.assertRaises(KeyError, get_mtz_flags, mtz )
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
