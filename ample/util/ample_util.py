'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

# Python modules
import cPickle
import logging
import os
import platform
import subprocess
import sys
import tarfile
import tempfile
import urllib
import zipfile

# our imports
import exit_util
import pdb_edit

from ample.constants import SHARE_DIR

CCP4_VERSION=None
SCRIPT_EXT = '.bat' if sys.platform.startswith('win') else '.sh'
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'

_logger = logging.getLogger(__name__)

# Reference string
references = """AMPLE: J. Bibby, R. M. Keegan, O. Mayans, M. D. Winn and D. J. Rigden.
AMPLE: a cluster-and-truncate approach to solve the crystal structures of small proteins using
rapidly computed ab initio models. (2012). Acta Cryst. D68, 1622-1631 [ doi:10.1107/S0907444912039194 ]

Routine phasing of coiled-coil protein crystal structures with AMPLE (2015). Thomas, J. M. H.,
Keegan, R. M., Bibby, J., Winn, M. D., Mayans, O. and Rigden, D. J. IUCrJ 2, 198-206.
[doi:10.1107/S2052252515002080] 

CCP4: Collaborative Computational Project, Number 4. (1994), The CCP4 Suite: Programs
for Protein Crystallography. Acta Cryst. D50, 760-763

MaxCluster: http://www.sbg.bio.ic.ac.uk/maxcluster/

MOLREP: A.A.Vagin & A.Teplyakov (1997) J. Appl. Cryst. 30, 1022-1025

MrBUMP: R.M.Keegan and M.D.Winn (2007) Acta Cryst. D63, 447-457

PHASER: McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D.,
Storoni, L.C. & Read, R.J. (2007)
Phaser crystallographic software J. Appl. Cryst. 40, 658-674

REFMAC: G.N. Murshudov, A.A.Vagin and E.J.Dodson, (1997) Refinement of Macromolecular
Structures by the Maximum-Likelihood Method. Acta Cryst. D53, 240-255

SCWRL: G. G. Krivov, M. V. Shapovalov, and R. L. Dunbrack, Jr. Improved prediction of protein
side-chain conformations with SCWRL4. Proteins (2009).

SHELXE: "Extending molecular-replacement solutions with SHELXE". Thorn, A. and Sheldrick, G. M. (2013),
Acta Crystallographica D, 69: 2251-2256. doi: 10.1107/S0907444913027534

SPICKER: Y. Zhang, J. Skolnick, SPICKER: Approach to clustering protein structures for
near-native model selection, Journal of Computational Chemistry, 2004 25: 865-871

Theseus: THESEUS: Maximum likelihood superpositioning and analysis of macromolecular structures.
Theobald, Douglas L. & Wuttke, Deborah S. (2006b) Bioinformatics 22(17):2171-2172 [Open Access]
Supplementary Materials for Theobald and Wuttke 2006b."""

# Header string
header ="""#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE - Ab Initio Modelling Molecular Replacement               #
#########################################################################

The authors of specific programs should be referenced where applicable:""" + \
"\n\n" + references + "\n\n"


survey_url = "http://goo.gl/forms/7xP9M4P81O"
footer = "\n" + """
#########################################################################
#***********************************************************************#
#*                          How did we do?                             *#
#*                                                                     *#
#*           Please follow this link and leave some feedback!          *#
#*                                                                     *#
#*                 {url}                      *#
#***********************************************************************#
#########################################################################
""".format(url=survey_url)


def ccp4_version():
    """Return the CCP4 version as a tuple"""
    global CCP4_VERSION
    if CCP4_VERSION is None:
        # Currently there seems no sensible way of doing this other then running a program and grepping the output
        cmd=['pdbcur']
        logf = tempfile.NamedTemporaryFile(delete=False)
        run_command(cmd, stdin="", logfile=logf.name)
        logf.seek(0) # rewind logfile
        tversion=None
        for i, line in enumerate(logf):
            if i > 20:break
            if line.startswith(' ### CCP4'):
                tversion=line.split()[2].rstrip(':')
                break
        
        logf.close()
        if not tversion: raise RuntimeError,"Cannot determine CCP4 version"
        vsplit = tversion.split('.')
        if len(vsplit) == 2:
            major = int(vsplit[0])
            minor =  int(vsplit[1])
            rev = '-1'
        elif len(vsplit) == 3:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = int(vsplit[2])
        else: raise RuntimeError,"Cannot split CCP4 version: {0}".format(tversion)
    os.unlink(logf.name)
    return (major,minor,rev)
    
def extract_models(amoptd, sequence=None, single=True, allsame=True):
    """Extract pdb files from a given tar/zip file or directory of pdbs"""
    
    filename = amoptd['models']
    directory = amoptd['models_dir']
    
    # If it's already a directory, just check it's valid   
    if os.path.isdir(filename):
        models_dir = filename
    else:
        # Here we are extracting from a file
        if not os.path.isfile(filename):
            msg="Cannot find models file: {0}".format(filename)
            exit_util.exit_error(msg)
            
        # we need a directory to extract into
        assert directory,"extractModels needs a directory path!"
        if not os.path.isdir(directory):
            os.mkdir(directory)
        models_dir = directory
        
        # See what sort of file this is:
        f,suffix=os.path.splitext(filename)
        if suffix in ['.gz','.bz']:
            f,s2=os.path.splitext(f)
            if s2 == '.tar': suffix=s2+suffix
        
        tsuffixes=['.tar.gz','.tgz','.tar.bz','.tbz']
        suffixes=tsuffixes + ['.zip']
        if suffix not in suffixes:
            msg="Do not know how to extract files from file: {0}\n Acceptable file types are: {1}".format(filename,suffixes)
            exit_util.exit_error(msg)
        if suffix in tsuffixes:
            files = extract_tar(filename, directory)
        else:
            files = extract_zip(filename, directory)
        
        # Assume anything with one member is quark decoys
        if len(files) == 1:
            quark_filename='alldecoy.pdb'
            f = os.path.basename(files[0])
            if not f == quark_filename:
                msg="Only found one member ({0}) in file: {1} and the name was not {2}\n".format(f, filename, quark_filename)
                msg+="If this file contains valid QUARK decoys, please email: ccp4@stfc.ac.uk"
                exit_util.exit_error(msg)
            # Now extract the quark pdb files from the monolithic file
            split_quark(files[0], models_dir)
            # We delete the quark_name file as otherwise we'll try and model it
            os.unlink(files[0])
            # If we've got quark models we don't want to modify the side chains as we only have polyalanine so we
            # set this here - horribly untidy as we should have one place to decide on side chains
            _logger.info('Found QUARK models in file: {0}'.format(filename))
            amoptd['quark_models'] = True
    
    if not pdb_edit.check_pdb_directory(models_dir, sequence=sequence, single=single, allsame=allsame):
        msg="Problem importing pdb files - please check the log for more information"
        exit_util.exit_error(msg)
    return models_dir

def extract_tar(filename, directory, suffixes=['.pdb']):
    # Extracting tarfile
    _logger.info('Extracting files from tarfile: {0}'.format(filename) )
    files = []
    with tarfile.open(filename,'r:*') as tf:
        memb = tf.getmembers()
        if not len(memb):
            msg='Empty archive: {0}'.format(filename)
            exit_util.exit_error(msg)
        for m in memb:
            if os.path.splitext(m.name)[1] in suffixes:
                # Hack to remove any paths
                m.name = os.path.basename(m.name)
                tf.extract(m, path=directory)
                files.append(os.path.join(directory, m.name))
    if not len(files):
        msg='Could not find any files with suffixes {0} in archive: {1}'.format(suffixes, filename)
        exit_util.exit_error(msg)
    return files

def extract_zip(filename, directory, suffixes=['.pdb']):
    # zip file extraction
    _logger.info('Extracting files from zipfile: {0}'.format(filename) )
    if not zipfile.is_zipfile(filename):
            msg='File is not a valid zip archive: {0}'.format(filename)
            exit_util.exit_error(msg)
    zipf = zipfile.ZipFile(filename)
    zif = zipf.infolist()
    if not len(zif):
        msg = 'Empty zip file: {0}'.format(filename)
        exit_util.exit_error(msg)
    files = []
    for f in zif:
        if os.path.splitext(f.filename)[1] in suffixes:
            # Hack to rewrite name 
            f.filename = os.path.basename(f.filename)
            zipf.extract(f, path=directory)
            files.append(os.path.join(directory, f.filename))
    if not len(files):
        msg = 'Could not find any files with suffixes {0} in zipfile: {1}'.format(suffixes,filename)
        exit_util.exit_error(msg)
    return files

def find_exe(executable, dirs=None):
    """Find the executable exename.
    Args:
    executable: the name of the program or the path to an existing executable
    dirs - additional directories to search for the location
    """
    _logger.debug('Looking for executable: {0}'.format(executable) )
    
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
        _logger.debug('Checking paths: {0}'.format(paths))
        
        for path in paths:
            exe_file = os.path.abspath(os.path.join(path, executable))   
            if is_exe(exe_file):
                _logger.debug( 'Found executable {0} in directory {1}'.format(executable,path) )
                found=True
                break
    
    if not found:
        raise Exception("Cannot find executable: {0}".format(executable))
    
    _logger.debug('find_exe found executable: {0}'.format(exe_file) )
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
        rcdir = amoptd['rcdir']
        _logger.info("Cannot find maxcluster binary in path so attempting to download it directory: {0}".format( rcdir )  )
        if not os.path.isdir( rcdir ):
            _logger.info("No ample rcdir found so creating in: {0}".format( rcdir ) )
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
                msg="Unrecognised system type: {0} {1}".format(sys.platform,bit)
                exit_util.exit_error(msg)
        elif sys.platform.startswith("darwin"):
            url = 'http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster_i686_32bit.bin'
            #OSX PPC: http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster_PPC_32bit.bin
        elif sys.platform.startswith("win"):
            url = 'http://www.sbg.bio.ic.ac.uk/~maxcluster/maxcluster.exe'
            maxcluster_exe = os.path.join( rcdir, 'maxcluster.exe' )
        else:
            msg="Unrecognised system type: {0}".format( sys.platform )
            exit_util.exit_error(msg)
        _logger.info("Attempting to download maxcluster binary from: {0}".format( url ) )
        try:
            urllib.urlretrieve( url, maxcluster_exe )
        except Exception, e:
            msg="Error downloading maxcluster executable: {0}\n{1}".format(url,e)
            exit_util.exit_error(msg)

        # make executable
        os.chmod(maxcluster_exe, 0o777)

    return maxcluster_exe

def ideal_helices(nresidues):
    ""
    include_dir = os.path.join(SHARE_DIR, 'include')
    names = [ 'polyala5', 'polyala10', 'polyala15', 'polyala20', 'polyala25',
              'polyala30', 'polyala35', 'polyala40' ]
    polya_lengths = [5,10,15,20,25,30,35,40]
    
    ensemble_options = {}
    ensembles_data = []
    pdbs = []
    for name, nres in zip(names, polya_lengths):
        ncopies = nresidues / nres
        if ncopies < 1: ncopies = 1
        ensemble_options[ name ] = { 'ncopies' : ncopies }
        pdb = os.path.join(include_dir,"{0}.pdb".format(name))
        # Needed for pyrvapi results
        ensembles_data.append( { 'name' : name,
                                'ensemble_pdb' : pdb,
                                'num_residues' : nres,
                                 } )
        pdbs.append(pdb)
        
    return pdbs, ensemble_options, ensembles_data

def is_exe(fpath):
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)

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

def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False, **kwargs):
    """Execute a command and return the exit code.

    We take care of outputting stuff to the logs and opening/closing logfiles

    Args:
    cmd - command to run as a list
    stdin - a string to use as stdin for the command
    logfile (optional) - the path to the logfile
    directory (optional) - the directory to run the job in (cwd assumed)
    dolog: bool - whether to output info to the system log
    """
    assert type(cmd) is list, "run_command needs a list!"
    if check:
        if not is_exe(cmd[0]): raise RuntimeError,"run_command cannot find executable: {0}".format(cmd[0])

    if not directory:  directory = os.getcwd()
    if dolog: _logger.debug("In directory {0}\nRunning command: {1}".format(directory, " ".join(cmd)))
    file_handle=False
    if logfile:
        if type(logfile)==file:
            file_handle=True
            logf=logfile
            logfile=os.path.abspath(logf.name)
        else:
            logfile = os.path.abspath(logfile)
            logf = open(logfile, "w")
        if dolog: _logger.debug("Logfile is: {0}".format(logfile))
    else:
        logf = tempfile.TemporaryFile()
        
    if stdin != None:
        stdinstr = stdin
        stdin = subprocess.PIPE

    # Windows needs some special treatment
    if os.name == "nt":
        kwargs.update( { 'bufsize': 0, 'shell' : "False" } )
    p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs)

    if stdin != None:
        p.stdin.write( stdinstr )
        p.stdin.close()
        if dolog: _logger.debug("stdin for cmd was: {0}".format( stdinstr ) )

    p.wait()
    if not file_handle: logf.close()
    return p.returncode

def saveAmoptd(amoptd):
    # Save results
    with open( amoptd['results_path'], 'w' ) as f:
        cPickle.dump( amoptd, f )
        _logger.info("Saved state as file: {0}\n".format( amoptd['results_path'] ) )
    return

def split_quark(dfile,directory):
    _logger.info("Extracting QUARK decoys from: {0} into {1}".format(dfile,directory))
    smodels = []
    with open(dfile,'r') as f:
        m=[]
        for line in f:
            if line.startswith("ENDMDL"):
                m.append(line)
                smodels.append(m)
                m = []
            else:
                m.append(line)
    if not len(smodels): raise RuntimeError,"Could not extract any models from: {0}".format(dfile)
    quark_models = []
    for i,m in enumerate(smodels):
        fpath = os.path.join(directory,"quark_{0}.pdb".format(i))
        with open(fpath,'w') as f:
            for l in m:
                # Need to reconstruct something sensible as from the coordinates on it's all quark-specific
                if l.startswith("ATOM"):
                    l = l[:54]+"  1.00  0.00              \n"
                f.write(l)
            quark_models.append(fpath)
            _logger.debug("Wrote: {0}".format(fpath))
        
    return quark_models

def tmpFileName():
    """Return a filename for a temporary file"""

    # Get temporary filenames
    t = tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True)
    tmp1 = t.name
    t.close()
    return tmp1
        

