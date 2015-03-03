'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# Python modules
import cPickle
import glob
import logging
import os
import platform
import re
import subprocess
import sys
import tarfile
import tempfile
import urllib
import zipfile

# External imports
import iotbx.pdb

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

def check_pdbs(directory,single=True):
    logger = logging.getLogger()
    if not os.path.isdir(directory):
        return False
    models=glob.glob(os.path.join(directory,"*.pdb"))
    if not len(models):
        return False
    
    if not single: return True
    
    # Check all models have 1 model and 1 chain
    def is_single(pdb):
        h=iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        return h.models_size()==1 and h.models()[0].chains_size()==1
    
    try:
        not_single = [ pdb for pdb in models if not is_single(pdb) ]
    except Exception,e:
        msg="Error processing pdbs in directory: {0}\n{1}".format(directory,e)
        logger.critical(msg)
        return False
    
    if len(not_single):
        msg="check_pdbs - the following pdb files have more than 1 chain/model:\n\n"
        for pdb in not_single: msg+="{0}\n".format(pdb)
        logger.critical(msg)
        return False
    else:
        logger.info("check_pdbs - pdb files all seem valid")
        return True

def extract_models(filename,directory=None):
    """Extract pdb files from a given tar/zip file or directory of pdbs"""
    
    logger = logging.getLogger()
    # If it's already a directory, just check it's valid   
    if os.path.isdir(filename):
        if not check_pdbs(filename):
            msg="Cannot extract pdb files from directory: {0}".format(filename)
            logger.critical(msg)
            raise RuntimeError,msg
        return filename

    # Here we are extracting from a file
    if not os.path.isfile(filename):
        msg="Cannot find models file: {0}".format(filename)
        logger.critical(msg)
        raise RuntimeError,msg
        
    # we need a directory to extract into
    assert directory,"extractModels needs a directory path!"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    models_dir=directory
    
    # See what sort of file this is:
    f,suffix=os.path.splitext(filename)
    if suffix in ['.gz','.bz']:
        f,s2=os.path.splitext(f)
        if s2 == '.tar': suffix=s2+suffix
    
    tsuffixes=['.tar.gz','.tgz','.tar.bz','.tbz']
    suffixes=tsuffixes + ['.zip']
    if suffix not in suffixes:
        msg="Do not know how to extract files from file: {0}\n Acceptable file types are: {1}".format(filename,suffixes)
        logger.critical(msg)
        raise RuntimeError,msg
    
    if suffix in tsuffixes:
        extract_tar(filename, directory)
    else:
        extract_zip(filename, directory)
        
    if not check_pdbs(models_dir):
        msg="Problem importing pdb files - please check the log for more information"
        logger.critical(msg)
        raise RuntimeError,msg
        
    return models_dir

def _extract_quark(tarfile,member,filename,models_dir):
    logger = logging.getLogger()
    # This is only acceptable if it is the quark decoys
    quark_name='alldecoy.pdb'
    if not member.name==quark_name:
        msg="Only found one member ({0}) in file: {1} and the name was not {2}\n".format(member.name,filename,quark_name)
        msg+="If this file contains valid QUARK decoys, please email: ccp4@stfc.ac.uk"
        logger.critical(msg)
        raise RuntimeError,msg
    
    # extract into current (work) directory
    tarfile.extract(member)
    
    # Now extract the quark pdb files from the monolithic file
    split_quark(member.name, models_dir)
    return

def extract_tar(filename,models_dir):
    # Extracting tarfile
    logger = logging.getLogger()
    logger.info('Extracting models from tarfile: {0}'.format(filename) )
    with tarfile.open(filename,'r:*') as tf:
        memb = tf.getmembers()
        if not len(memb):
            msg='Empty archive: {0}'.format(filename)
            logger.critical(msg)
            raise RuntimeError,msg
        if len(memb) == 1:
            # Assume anything with one member is quark decoys
            logger.info('Checking if file contains quark decoys'.format(filename))
            _extract_quark(tf,memb[0],filename,models_dir)
        else:
            got=False
            for m in memb:
                if os.path.splitext(m.name)[1] == '.pdb':
                    # Hack to remove any paths
                    m.name=os.path.basename(m.name)
                    tf.extract(m,path=models_dir)
                    got=True
            if not got:
                msg='Could not find any pdb files in archive: {0}'.format(filename)
                logger.critical(msg)
                raise RuntimeError,msg
    return

def extract_zip(filename,models_dir,suffix='.pdb'):
    # zip file extraction
    logger = logging.getLogger()
    logger.info('Extracting models from zipfile: {0}'.format(filename) )
    if not zipfile.is_zipfile(filename):
            msg='File is not a valid zip archive: {0}'.format(filename)
            logger.critical(msg)
            raise RuntimeError,msg
    zipf=zipfile.ZipFile(filename)
    zif=zipf.infolist()
    if not len(zif):
        msg='Empty zip file: {0}'.format(filename)
        logger.critical(msg)
        raise RuntimeError,msg
    got=False
    for f in zif:
        if os.path.splitext(f.filename)[1] == suffix:
            # Hack to rewrite name 
            f.filename=os.path.basename(f.filename)
            zipf.extract(f, path=models_dir)
            got=True
    if not got:
        msg='Could not find any pdb files in zipfile: {0}'.format(filename)
        logger.critical(msg)
        raise RuntimeError,msg    
    return

def find_exe(executable, dirs=None):
    """Find the executable exename.

    Args:
    executable: the name of the program or the path to an existing executable
    dirs - additional directories to search for the location
    
    """

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

def filename_append(filename=None, astr=None,directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split( filename )
    name, suffix = os.path.splitext( fname )
    name  =  name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join( directory, name )

def find_maxcluster(amoptd):
    """Return path to maxcluster binary.
    If we can't find one in the path, we create a $HOME/.ample
    directory and downlod it to there
    """

    if amoptd['maxcluster_exe'] and is_exe(amoptd['maxcluster_exe']):
        return amoptd['maxcluster_exe']

    if not amoptd['maxcluster_exe']:
        if sys.platform.startswith("win"):
            amoptd['maxcluster_exe']='maxcluster.exe'
        else:
            amoptd['maxcluster_exe']='maxcluster'
    
    try:
        maxcluster_exe = find_exe(amoptd['maxcluster_exe'], dirs=[ amoptd['rcdir'] ] )
    except Exception:
        # Cannot find so we need to try and download it
        logger = logging.getLogger()
        rcdir = amoptd['rcdir']
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

def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

def make_workdir(work_dir, ccp4_jobid=None, rootname='AMPLE_'):
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

def saveAmoptd(amoptd):
    # Save results
    with open( amoptd['results_path'], 'w' ) as f:
        cPickle.dump( amoptd, f )
        logging.info("Saved results as file: {0}\n".format( amoptd['results_path'] ) )
    return

def split_quark(dfile,directory):
    logger = logging.getLogger()
    logger.info("Extracting QUARK decoys from: {0} into {1}".format(dfile,directory))
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
    if not len(smodels): raise RuntimeError,"Could not extract any models from: {0}".format(dfile)
    for i,m in enumerate(smodels):
        fpath=os.path.join(directory,"quark_{0}.pdb".format(i))
        with open(fpath,'w') as f:
            for l in m:
                # Need to reconstruct something sensible as from the coordinates on it's all quark-specific
                if l.startswith("ATOM"):
                    l=l[:54]+"  1.00  0.00              \n"
                f.write(l)
        logger.debug("Wrote: {0}".format(fpath))
    return

def tmpFileName():
    """Return a filename for a temporary file"""

    # Get temporary filenames
    t = tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True)
    tmp1 = t.name
    t.close()
    return tmp1
