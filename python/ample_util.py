'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# Python modules
import logging
import os
import platform
import re
import subprocess
import sys
import tarfile
import tempfile
import urllib
import unittest

# Our modules
import cif_parser

# MRBUMP modules
#sys.path.append(os.path.join(os.environ["CCP4"], "share", "mrbump", "include", "file_info")) # For MTZ_parse
sys.path.insert(0, "/Users/jmht/Documents/AMPLE/programs/mrbump/include/file_info")
import MTZ_parse

# Reference string
references = """AMPLE: J. Bibby, R. M. Keegan, O. Mayans, M. D. Winn and D. J. Rigden.
AMPLE: a cluster-and-truncate approach to solve the crystal structures of small proteins using
rapidly computed ab initio models. (2012). Acta Cryst. D68, 1622-1631 [ doi:10.1107/S0907444912039194 ]

CCP4: Collaborative Computational Project, Number 4. (1994), The CCP4 Suite: Programs
for Protein Crystallography. Acta Cryst. D50, 760-763

SHELX: "A short history of SHELX". Sheldrick, G.M. (2008). Acta Cryst. A64, 112-122

SCWRL: G. G. Krivov, M. V. Shapovalov, and R. L. Dunbrack, Jr. Improved prediction of protein
side-chain conformations with SCWRL4. Proteins (2009).

MaxCluster: http://www.sbg.bio.ic.ac.uk/maxcluster/

MOLREP: A.A.Vagin & A.Teplyakov (1997) J. Appl. Cryst. 30, 1022-1025

MrBUMP: R.M.Keegan and M.D.Winn (2007) Acta Cryst. D63, 447-457

PHASER: McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D.,
Storoni, L.C. & Read, R.J. (2007)
Phaser crystallographic software J. Appl. Cryst. 40, 658-674

REFMAC: G.N. Murshudov, A.A.Vagin and E.J.Dodson, (1997) Refinement of Macromolecular
Structures by the Maximum-Likelihood Method. Acta Cryst. D53, 240-255

SPICKER: Y. Zhang, J. Skolnick, SPICKER: Approach to clustering protein structures for
near-native model selection, Journal of Computational Chemistry, 2004 25: 865-871

Theseus: THESEUS: Maximum likelihood superpositioning and analysis of macromolecular structures.
Theobald, Douglas L. & Wuttke, Deborah S. (2006b) Bioinformatics 22(17):2171-2172 [Open Access]
Supplementary Materials for Theobald and Wuttke 2006b."""

# Header string
header ="""#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE -Ab Initio Modelling Molecular Replacement (Beta version) #
#########################################################################
The authors of specific programs should be referenced where applicable:""" + \
"\n\n" + references + "\n\n"

def extractFile(tarArchive,fileName,directory=None):
    """Extract a file from a tar.gz archive into the directory and return the name of the file"""
    with tarfile.open(tarArchive,'r:*') as tf:
        m = tf.getmembers()
        if not len(m):
            raise RuntimeError,'Empty archive: {0}'.format(tarArchive)
        if len(m)==1 and m[0].name==fileName:
            tf.extractall(path=directory)
            return m[0].name
        else:
            return False

def find_exe( executable, dirs=None ):
    """Find the executable exename.

    Args:
    executable: the name of the program or the path to an existing executable
    dirs - additional directories to search for the location
    
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    logger = logging.getLogger()
    logger.debug('Looking for executable: {0}'.format(executable) )
    
    exe_file=None
    found=False
    if is_exe(executable):
        exe_file=os.path.abspath(executable)
        found=True
    else:
        # If the user has given a path we just take the name part
        fpath,fname = os.path.split(executable)
        if fname:
            executable=fname
            
        # By default we search in the system PATH and add any additional user given paths here
        paths = os.environ["PATH"].split(os.pathsep)
        if dirs:
            paths += dirs
        logger.debug('Checking paths: {0}'.format(paths))
        for path in paths:
            exe_file = os.path.abspath(os.path.join(path, executable))
            if is_exe(exe_file):
                logger.debug( 'Found executable {0} in directory {1}'.format(executable,path) )
                found=True
                break
    
    if not found:
        raise Exception("Cannot find executable: {0}".format(executable))
    
    logger.debug('find_exe found executable: {0}'.format(exe_file) )
    return exe_file

def filename_append( filename=None, astr=None,directory=None, separator="_",  ):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split( filename )
    name, suffix = os.path.splitext( fname )
    name  =  name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join( directory, name )

def find_maxcluster( amopt ):
    """Return path to maxcluster binary.
    If we can't find one in the path, we create a $HOME/.ample
    directory and downlod it to there
    """

    if not amopt.d['maxcluster_exe']:
        if sys.platform.startswith("win"):
            amopt.d['maxcluster_exe']='maxcluster.exe'
        else:
            amopt.d['maxcluster_exe']='maxcluster'
    
    try:
        maxcluster_exe = find_exe(amopt.d['maxcluster_exe'], dirs=[ amopt.d['rcdir'] ] )
    except Exception:
        # Cannot find so we need to try and download it
        logger = logging.getLogger()
        rcdir = amopt.d['rcdir']
        logger.info("Cannot find maxcluster binary in path so attempting to download it directory: {0}".format( rcdir )  )
        if not os.path.isdir( rcdir ):
            logger.info("No ample rcdir found so creating in: {0}".format( rcdir ) )
            os.mkdir( rcdir )
        url = None
        maxcluster_exe = os.path.join( rcdir, 'maxcluster' )
        if sys.platform.startswith("linux"):
            bit=platform.architecture()[0]
            if bit=='64bit':
                url='http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster64bit'
            elif bit=='32bit':
                url='http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster'
            else:
                raise RuntimeError,"Unrecognised system type: {0} {1}".format(sys.platform,bit)
        elif sys.platform.startswith("darwin"):
            url = 'http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster_i686_32bit.bin'
            #OSX PPC: http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster_PPC_32bit.bin
        elif sys.platform.startswith("win"):
            url = 'http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster.exe'
            maxcluster_exe = os.path.join( rcdir, 'maxcluster.exe' )
        else:
            raise RuntimeError,"Unrecognised system type: {0}".format( sys.platform )

        logger.info("Attempting to download maxcluster binary from: {0}".format( url ) )
        try:
            urllib.urlretrieve( url, maxcluster_exe )
        except Exception, e:
            logger.critical("Error downloading maxcluster executable: {0}\n{1}".format( url, e ) )
            raise

        # make executable
        os.chmod(maxcluster_exe, 0o777)

    return maxcluster_exe

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

def processReflectionFile( amoptd ):
    """Make sure we have a valid mtz file. If necessary convert a given cif file.
       Set the mtz variable in the given amoptd to the reflection file to use
       Return True if it all worked or raise an exception if it failed
    """

    # We've been given a sf_cif so convert to mtz
    if amoptd['sf_cif']:
        if not os.path.isfile( amoptd['sf_cif'] ):
            logging.critical("Cannot find sf_cif file: {0}".format( amoptd['sf_cif'] ) )
            sys.exit(1)

        cp = cif_parser.CifParser()
        amoptd['mtz'] = cp.sfcif2mtz( amoptd['sf_cif'] )

    # Now have an mtz so check it's valid
    if not amoptd['mtz'] or not os.path.isfile( amoptd['mtz'] ):
        logging.critical("Cannot find MTZ file: {0}".format( amoptd['mtz'] ) )
        sys.exit(1)

    # Run mtzdmp to get file info
    mtzp = MTZ_parse.MTZ_parse()
    mtzp.run_mtzdmp( amoptd['mtz'] )

    # If flags given check they are present in the file
    for flag in ['F', 'SIGF', 'FREE']:
        if amoptd[flag] and amoptd[flag] not in mtzp.col_labels:
            logging.critical("Cannot find given {0} label {1} in mtz file {2}".format( flag, amoptd[flag], amoptd['mtz'] ) )
            sys.exit(1)

    # If any of the flags aren't given we get them from the file
    if not amoptd['F']:
        amoptd['F']  = mtzp.F
    if not amoptd['SIGF']:
        amoptd['SIGF']  = mtzp.SIGF
    if not amoptd['FREE']:
        amoptd['FREE']  = mtzp.FreeR_flag

    # Make sure we've found something
    for flag in ['F', 'SIGF' ]:
        if amoptd[flag] is None or not len(amoptd[flag]) or amoptd[flag] not in mtzp.col_labels:
            logging.critical("Cannot find any {0} label in mtz file {1}".format( flag, amoptd['mtz'] ) )
            sys.exit(1)

    # All flags ok so just check the RFREE flag is valid
    if not mtzp.checkRFREE(FreeR_flag=amoptd['FREE']):
        # If not run uniqueify
        logging.warning("Cannot find a valid FREE flag - running uniquefy to generate column with RFREE data." )
        amoptd['mtz'] = uniqueify( amoptd['mtz'], directory=amoptd['work_dir'] )

        # Check file and get new FREE flag
        mtzp.run_mtzdmp( amoptd['mtz'] )
        amoptd['FREE']  = mtzp.FreeR_flag

        # hopefully unnecessary check
        assert mtzp.checkRFREE(FreeR_flag=amoptd['FREE'])

    return True

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

    # Windows needs some special treatment
    kwargs = {}
    if os.name == "nt":
        kwargs = { 'bufsize': 0, 'shell' : "False" }
    p = subprocess.Popen( cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs )

    if stdin != None:
        p.stdin.write( stdinstr )
        p.stdin.close()
        if dolog:
            logging.debug("stdin for cmd was: {0}".format( stdinstr ) )

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

def splitQuark(dfile,directory='quark_models'):
    smodels=[]
    with open(dfile,'r') as f:
        m=[]
        for line in f:
            if line.startswith("ENDMDL"):
                m.append(line)
                smodels.append(m)
                m=[]
            else:
                m.append(line)

    if not len(smodels):
        raise RuntimeError,"Could not extract any models from: {0}".format(dfile)

    if not os.path.isdir(directory):
        os.mkdir(directory)

    for i,m in enumerate(smodels):
        fpath=os.path.join(directory,"quark_{0}.pdb".format(i))
        with open(fpath,'w') as f:
            for l in m:
                # Need to reconstruct something sensible as from the coordinates on it's all quark-specific
                if l.startswith("ATOM"):
                    l=l[:54]+"  1.00  0.00              \n"
                f.write(l)
        logging.debug("Wrote: {0}".format(fpath))

    return

def tmpFileName():
    """Return a filename for a temporary file"""

    # Get temporary filenames
    t = tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True)
    tmp1 = t.name
    t.close()
    return tmp1

def uniqueify(mtzPath,directory=None):
    """Run uniqueify on mtz file to generate RFREE data column"""

    mtzUnique = filename_append(mtzPath, "uniqueify", directory=directory)

    cmd = ['uniqueify', mtzPath, mtzUnique]

    #uniqueify {-p fraction} mydata.mtz.
    logfile = os.path.join( os.getcwd(), "uniqueify.log" )
    retcode = run_command(cmd, logfile=logfile)
    if retcode != 0:
        msg = "Error running sfcif2mtz. Check the logfile: {0}".format(logfile)
        logging.critical(msg)
        raise RuntimeError, msg

    return mtzUnique

class TestUtil(unittest.TestCase):
    """
    Unit test
    """

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ampleDir = os.sep.join( paths[ : -1 ] )
        self.testfilesDir = os.sep.join( paths[ : -1 ] + [ 'tests', 'testfiles' ] )

        return

    def testProcessReflectionFile(self):
        """Get MTZ flags"""


        mtz = os.path.join( self.ampleDir, "examples", "toxd-example" , "1dtx.mtz" )


        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
             }

        processReflectionFile( d )

        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        return

    def testProcessReflectionFileNORFREE(self):
        """Get MTZ flags"""

        mtz = os.path.join( self.testfilesDir, "2uui_sigmaa.mtz" )

        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
             }

        processReflectionFile( d )

        self.assertEqual( 'F', d['F'], "Correct F")
        self.assertEqual( 'SIGF', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        return

    def testProcessReflectionFileCIF(self):
        """Get MTZ flags"""

        cif = os.path.join( self.testfilesDir, "1x79-sf.cif" )

        d = { 'mtz'    : None,
              'sf_cif' : cif,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
             }

        processReflectionFile( d )

        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        return
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
