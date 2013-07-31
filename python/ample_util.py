'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# Python modules
import logging
import os
import re
import subprocess
import sys
import tempfile
import unittest

# Our modules
import MTZParse

# Reference string
references = """AMPLE: To be added

SHELX: is used: "A short history of SHELX". Sheldrick, G.M. (2008). Acta Cryst. A64, 112-122

SCWRL: G. G. Krivov, M. V. Shapovalov, and R. L. Dunbrack, Jr. Improved prediction of protein
side-chain conformations with SCWRL4. Proteins (2009).

Theseus: THESEUS: Maximum likelihood superpositioning and analysis of macromolecular structures.
Theobald, Douglas L. & Wuttke, Deborah S. (2006b) Bioinformatics 22(17):2171-2172 [Open Access]
Supplementary Materials for Theobald and Wuttke 2006b.

MrBUMP: R.M.Keegan and M.D.Winn (2007) Acta Cryst. D63, 447-457

CCP4: Collaborative Computational Project, Number 4. (1994), The CCP4 Suite: Programs
for Protein Crystallography. Acta Cryst. D50, 760-763

MOLREP: A.A.Vagin & A.Teplyakov (1997) J. Appl. Cryst. 30, 1022-1025

PHASER: McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D.,
Storoni, L.C. & Read, R.J. (2007)
Phaser crystallographic software J. Appl. Cryst. 40, 658-674

REFMAC: G.N. Murshudov, A.A.Vagin and E.J.Dodson, (1997) Refinement of Macromolecular
Structures by the Maximum-Likelihood Method. Acta Cryst. D53, 240-255

SPICKER: Y. Zhang, J. Skolnick, SPICKER: Approach to clustering protein structures for
near-native model selection, 
Journal of Computational Chemistry, 2004 25: 865-871

MaxCluster: http://www.sbg.bio.ic.ac.uk/maxcluster/"""

# Header string
header ="""#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE -Ab Initio Modelling Molecular Replacement (Beta version) #
#########################################################################
The authors of specific programs should be referenced where applicable:""" + \
"\n\n" + references + "\n\n"

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
    logger.debug('Looking for executable: {0}'.format(exename) )
    
    if varname:
        logger.debug( 'Using executable path from command-line {0}'.format(varname) )
        exepath = which(varname)
        if not exepath:
            msg = "Cannot find valid {0} executable from given path: {1}".format(exename,varname)
            logger.critical(msg)
            raise RuntimeError,msg
    else:
        logger.debug( 'No path for {0} given on the command line, looking in PATH'.format(exename) )
        exepath = which(exename)
        if not exepath:
            msg = "Cannot find executable {0} in PATH. Please give the path to {0}".format(exename)
            logger.critical(msg)
            raise RuntimeError,msg

    logger.debug( "Using executable {0}".format( exepath ) )
    return exepath


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
    
    try:
        F =  md.colLabels['F']
    except KeyError:
        raise RuntimeError,"Cannot find column label: {0} in MTZ file".format('F')
    try:
        SIGF =  md.colLabels['SIGF']
    except KeyError:
        raise RuntimeError,"Cannot find column label: {0} in MTZ file".format('SIGF')
    try:
        FREE =  md.colLabels['FREE']
    except KeyError:
        raise RuntimeError,"Cannot find column label: {0} in MTZ file".format('FREE')
    
    return (F, SIGF, FREE)
    
##End get_mtz_flags

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

########

def make_workdir(work_dir, ccp4_jobid=None, rootname='ROSETTA_MR_'):
    """
    Make a work directory rooted at work_dir and return its path
    """
    
    if ccp4_jobid:
        dname = os.path.join( work_dir, rootname + str(ccp4_jobid) )
        if os.path.exists(dname):
            raise RuntimeError,"There is an existing AMPLE CCP4 work directory: {0}\nPlease delete/move it aside."
        os.mkdir(dname)
        return dname
    
    run_inc = 0
    run_making_done = False
    while not run_making_done:
        if not os.path.exists(work_dir + os.sep + rootname + str(run_inc)):
            run_making_done = True
            os.mkdir(work_dir + os.sep +rootname + str(run_inc))
        run_inc += 1
    work_dir = work_dir + os.sep + rootname + str(run_inc - 1)
    return work_dir

def run_command( cmd, logfile=None, directory=None, dolog=True, stdin=None ):
    """Execute a command and return the exit code.
    
    We take care of outputting stuff to the logs and opening/closing logfiles
    
    Args:
    cmd - command to run as a list
    stdin - a string to use as stdin for the command
    logfile (optional) - the path to the logfile
    directory (optional) - the directory to run the job in (cwd assumed)
    dolog: bool - whether to output info to the system log
    """
    
    assert type(cmd) is list
    
    if not directory:
        directory = os.getcwd()
        
    if dolog:
        logging.debug("In directory {0}\nRunning command: {1}".format( directory, " ".join(cmd)  ) )
    
    if logfile:
        if dolog:
            logging.debug("Logfile is: {0}".format( logfile ) )
        logf = open( logfile, "w" )
    else:
        logf = tempfile.TemporaryFile()
        
    if stdin != None:
        stdinstr = stdin
        stdin = subprocess.PIPE
    
    p = subprocess.Popen( cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory )
    
    if stdin != None:
        p.stdin.write( stdinstr )
        p.stdin.close()
    
    p.wait()
    logf.close()
    
    return p.returncode

def setup_logging():
    """
    Set up the various log files/console logging
    and return the logger
    """
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # create file handler and set level to debug
    fl = logging.FileHandler("debug.log")
    fl.setLevel(logging.DEBUG)
    
    # create formatter for fl
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # add formatter to fl
    fl.setFormatter(formatter)
    
    # add fl to logger
    logger.addHandler(fl)
    
    # Now create console logger for outputting stuff
    # create file handler and set level to debug
    # Seems they changed the api in python 2.6->2.7
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(strm=sys.stdout)

    cl.setLevel(logging.INFO)
    
    # create formatter for fl
    # Always add a blank line after every print
    formatter = logging.Formatter('%(message)s\n')
    
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
