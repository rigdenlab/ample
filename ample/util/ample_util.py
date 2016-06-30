'''
Various miscellaneous functions.
Might end up somewhere else at somepoint.
'''

import cPickle
import logging
import os
import platform
import subprocess
import sys
import tarfile
import tempfile
import urllib
import warnings
import zipfile

# our imports
import exit_util
import pdb_edit

from ample.constants import SHARE_DIR

CCP4_VERSION=None
SCRIPT_EXT = '.bat' if sys.platform.startswith('win') else '.sh'
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'

LOGGER = logging.getLogger(__name__)


def ccp4_version():
    """
    Get CCP4 version as a tuple

    Returns
    -------
    version : tuple
       Major, minor, and revision number
    """
    global CCP4_VERSION
    if CCP4_VERSION is None:
        # Currently there seems no sensible way of doing this other then running a program and grepping the output
        pdbcur = 'pdbcur' + EXE_EXT
        log_fname = tmp_file_name(delete=False)
        run_command([pdbcur], stdin="", logfile=log_fname)
        tversion = None

        with open(log_fname, 'r') as logfh:
            for i, line in enumerate(logfh.readlines()):
                if i > 20:
                    break
                if line.startswith(' ### CCP4'):
                    tversion = line.split()[2].rstrip(':')
                    break

        if not tversion:
            raise RuntimeError("Cannot determine CCP4 version")
        vsplit = tversion.split('.')
        if len(vsplit) == 2:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = '-1'
        elif len(vsplit) == 3:
            major = int(vsplit[0])
            minor = int(vsplit[1])
            rev = int(vsplit[2])
        else:
            raise RuntimeError("Cannot split CCP4 version: {0}".format(tversion))
    os.unlink(log_fname)
    return (major,minor,rev)

def construct_references():
    """
    Construct the reference list

    Description
    -----------
    This is somewhat a very basic BibTex parser. It is under no circumstances foolproof but rather
    the bare minimum for what is required in AMPLE.

    If a more sophisticated parser is required, refer to one of the python packages available online.

    Returns
    -------
    references : str
    """

    # Get the filename and check we can use it
    ref_fname = os.path.join(SHARE_DIR, "include", "references.bib")
    if not is_file(ref_fname):
        msg = "Cannot find BibTex file containing references. " \
              "Please determine them yourself and cite AMPLE."
        return msg

    articles = []
    article = {}
    entry = False

    with open(ref_fname, "r") as fhin:
        for line in fhin.readlines():

            # Make sure line is not empty
            line = line.strip()
            if not line:
                continue

            # Beginning of all BibTex entries
            if line.startswith("@"):
                entry = True
                article = {}                # Reset the article dictionary
            # End of of all BibTex entries
            elif line.startswith("}") and line.endswith("}"):
                entry = False
                articles.append(article)    # Save the article dictionary

            if entry and not line.startswith("@"):

                line = line.replace("{", "")
                line = line.replace("}", "")
                line = line.split("=")

                line = [line[0].strip(), line[1].strip().rstrip(",").replace("\"", "")]

                # Do some formatting of the data so we can have a properly formatted string
                # We might also want to - at some point ever - use this data otherhow so we
                # this allows us to do it
                if line[0].lower() == "author":
                    tmp = line[1].split(" and ")[0]
                    line[1] = tmp.split(",")[0] + " et al."
                if line[0].lower() == "volume":
                    line[1] = int(line[1])
                if line[0].lower() == "year":
                    line[1] = int(line[1])
                if line[0].lower() == "number":
                    line[1] = int(line[1])
                if line[0].lower() == "pages":
                    line[1] = line[1].replace("--", "-")

                # Save the data to the article entry
                article[line[0]] = line[1]

    # Somewhat a template of how we want to display each article
    template = "* {label}: {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
    references = [template.format(**article) for article in articles]

    return (os.linesep*2).join(references)

def extract_models(amoptd, sequence=None, single=True, allsame=True):
    """Check a directory of pdbs or extract pdb files from a given tar/zip file or directory of pdbs
    and set the amoptd['models_dir'] entry with the directory of unpacked/validated pdbs
    """
    
    filename = amoptd['models']
    models_dir = amoptd['models_dir']
    
    # If it's already a models_dir, just check it's valid   
    if os.path.isdir(filename):
        models_dir = filename
    else:
        # Here we are extracting from a file
        if not os.path.isfile(filename):
            msg="Cannot find models file: {0}".format(filename)
            exit_util.exit_error(msg)
            
        # we need a models_dir to extract into
        assert models_dir,"extractModels needs a models_dir path!"
        if not os.path.isdir(models_dir):
            os.mkdir(models_dir)
        models_dir = models_dir
        
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
            files = extract_tar(filename, models_dir)
        else:
            files = extract_zip(filename, models_dir)
        
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
            LOGGER.info('Found QUARK models in file: {0}'.format(filename))
            amoptd['quark_models'] = True
    
    if not pdb_edit.check_pdb_directory(models_dir, sequence=sequence, single=single, allsame=allsame):
        msg="Problem importing pdb files - please check the log for more information"
        exit_util.exit_error(msg)
    
    amoptd['models_dir'] = models_dir
    return

def extract_tar(filename, directory, suffixes=['.pdb']):
    # Extracting tarfile
    LOGGER.info('Extracting files from tarfile: {0}'.format(filename) )
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
    LOGGER.info('Extracting files from zipfile: {0}'.format(filename) )
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
    LOGGER.debug('Looking for executable: {0}'.format(executable) )
    
    exe_file=None
    found=False
    if is_exe(executable):
        exe_file=os.path.abspath(executable)
        found=True
    else:
        # If the user has given a path we just take the name part
        _, fname = os.path.split(executable)
        if fname:
            executable=fname
            
        # By default we search in the system PATH and add any additional user given paths here
        paths = os.environ["PATH"].split(os.pathsep)
        if dirs:
            paths += dirs
        LOGGER.debug('Checking paths: {0}'.format(paths))
        
        for path in paths:
            exe_file = os.path.abspath(os.path.join(path, executable))   
            if is_exe(exe_file):
                LOGGER.debug( 'Found executable {0} in directory {1}'.format(executable,path) )
                found=True
                break
    
    if not found:
        raise Exception("Cannot find executable: {0}".format(executable))
    
    LOGGER.debug('find_exe found executable: {0}'.format(exe_file) )
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
        LOGGER.info("Cannot find maxcluster binary in path so attempting to download it directory: {0}".format( rcdir )  )
        if not os.path.isdir( rcdir ):
            LOGGER.info("No ample rcdir found so creating in: {0}".format( rcdir ) )
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
        LOGGER.info("Attempting to download maxcluster binary from: {0}".format( url ) )
        try:
            urllib.urlretrieve( url, maxcluster_exe )
        except Exception, e:
            msg="Error downloading maxcluster executable: {0}\n{1}".format(url,e)
            exit_util.exit_error(msg)

        # make executable
        os.chmod(maxcluster_exe, 0o777)

    return maxcluster_exe

def ideal_helices(nresidues):
    """
    Get some ideal helices

    Parameters
    ----------
    nresidues : int
       Number of residues to be used

    Returns
    -------
    pdbs : list
    ensemble_options : dict
    ensembles_data : list
    """
    include_dir = os.path.join(SHARE_DIR, 'include')
    names = ['polyala5', 'polyala10', 'polyala15', 'polyala20', 'polyala25',
             'polyala30', 'polyala35', 'polyala40']
    polya_lengths = [5, 10, 15, 20, 25, 30, 35, 40]
    
    ensemble_options = {}
    ensembles_data = []
    pdbs = []
    for name, nres in zip(names, polya_lengths):
        ncopies = nresidues / nres
        if ncopies < 1: ncopies = 1
        ensemble_options[name] = {'ncopies' : ncopies}
        pdb = os.path.join(include_dir, "{0}.pdb".format(name))
        # Needed for pyrvapi results
        ensembles_data.append({'name': name,
                               'ensemble_pdb': pdb,
                               'num_residues': nres,
                                })
        pdbs.append(pdb)
        
    return pdbs, ensemble_options, ensembles_data

def is_exe(fpath):
    """
    Check if an executable exists

    Parameters
    ----------
    fpath : str
       The path to the executable
    """
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)

def is_file(fpath):
    """
    Check if a file exists

    Parameters
    ----------
    fpath : str
       The path to the file
    """
    return fpath and os.path.isfile(fpath) and \
           os.access(fpath, os.R_OK) and os.stat(fpath).st_size > 0

def make_workdir(work_dir, ccp4_jobid=None, rootname='AMPLE_'):
    """
    Make a work directory rooted at work_dir and return its path

    Parameters
    ----------
    work_dir : str
       The path to a working directory
    ccp4_jobid : int
       CCP4-assigned job identifier
    rootname : str
        Base name of the AMPLE directory

    Returns
    -------
    work_dir : str
       The path to the working directory
    """

    if ccp4_jobid:
        dname = os.path.join(work_dir, rootname + str(ccp4_jobid))
        if os.path.exists(dname):
            raise RuntimeError("There is an existing AMPLE CCP4 work directory: {0}\n"
                               "Please delete/move it aside.")
        os.mkdir(dname)
        return dname

    run_inc = 0
    run_making_done = False
    while not run_making_done:
        if not os.path.exists(work_dir + os.sep + rootname + str(run_inc)):
            run_making_done = True
            os.mkdir(work_dir + os.sep + rootname + str(run_inc))
        run_inc += 1
    work_dir = work_dir + os.sep + rootname + str(run_inc - 1)
    return work_dir

def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False, **kwargs):
    """
    Execute a command and return the exit code.

    Parameters
    ----------
    cmd : list
       Command to run as a list
    stdin : str
       Stdin for the command
    logfile : str
       The path to the logfile
    directory : str
       The directory to run the job in (cwd assumed)
    dolog : bool
       Whether to output info to the system log

    Returns
    -------
    returncode : int
       Subprocess exit code

    Notes
    -----
    We take care of outputting stuff to the logs and opening/closing logfiles
    """
    assert type(cmd) is list, "run_command needs a list!"
    if check and not is_exe(cmd[0]):
        raise RuntimeError("run_command cannot find executable: {0}".format(cmd[0]))

    if not directory:
        directory = os.getcwd()

    if dolog:
        LOGGER.debug("In directory {0}".format(directory))
        LOGGER.debug("Running command: {0}".format(" ".join(cmd)))
        if kwargs:
            LOGGER.debug("kwargs are: {0}".format(kwargs))

    file_handle = False
    if logfile:
        if type(logfile) == file:
            file_handle = True
            logf = logfile
            logfile = os.path.abspath(logf.name)
        else:
            logfile = os.path.abspath(logfile)
            logf = open(logfile, "w")
        if dolog: LOGGER.debug("Logfile is: {0}".format(logfile))
    else:
        logf = tmp_file_name()
        
    if stdin != None:
        stdinstr = stdin
        stdin = subprocess.PIPE

    # Windows needs some special treatment
    if os.name == "nt":
        kwargs.update({'bufsize': 0, 'shell' : "False"})
    p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs)

    if stdin != None:
        p.stdin.write(stdinstr)
        p.stdin.close()
        if dolog: LOGGER.debug("stdin for cmd was: {0}".format(stdinstr))

    p.wait()
    if not file_handle:
        logf.close()

    return p.returncode

def read_amoptd(amoptd_fname):
    """
    Read a PICKLE-formatted AMPLE options file

    Parameters
    ----------
    amoptd_fname : str
       The path to the PICKLE-formatted AMPLE options file

    Returns
    -------
    amoptd : dict
       AMPLE options from saved state
    """
    if not is_file(amoptd_fname):
        raise RuntimeError("Something is wrong with your AMPLE options "
                           "file: {0}\n".format(amoptd_fname))

    with open(amoptd_fname, 'r') as f:
        amoptd = cPickle.load(f)
        LOGGER.info("Loaded state from file: {0}\n".format(amoptd['results_path']))
    return amoptd

def saveAmoptd(*args):
    """
    Save AMPLE options to a PICKLE-formatted file

    See Also
    --------
    save_amoptd

    Warnings
    --------
    This function was deprecated and will be removed in future releases. Please use ``save_amoptd()`` instead.
    """
    msg = "This function was deprecated and will be removed in future release"
    warnings.warn(msg, DeprecationWarning, stacklevel=2)
    save_amoptd(*args)
    return

def save_amoptd(amoptd):
    """
    Save AMPLE options to a PICKLE-formatted file

    Parameters
    ----------
    amoptd : dict
       AMPLE options from saved state
    """
    # Save results
    with open(amoptd['results_path'], 'w') as f:
        cPickle.dump(amoptd, f)
        LOGGER.info("Saved state as file: {0}\n".format(amoptd['results_path']))
    return

def split_quark(*args):
    """
    Split a single PDB with multiple models in individual PDB files

    See Also
    --------
    split_models
    """
    return split_models(*args)

def split_models(dfile, directory):
    """
    Split a single PDB with multiple models in individual PDB files

    Parameters
    ----------
    dfile : str
       Single PDB file with multiple model entries
    directory : str
       Directory to extract the PDB files to

    Returns
    -------
    extracted_models : list
       List of PDB files for all models

    TODO
    ----
    * Use the CCTBX library to perform this step
    """
    LOGGER.info("Extracting decoys from: {0} into {1}".format(dfile, directory))
    smodels = []
    with open(dfile, 'r') as f:
        m = []
        for line in f:
            if line.startswith("ENDMDL"):
                m.append(line)
                smodels.append(m)
                m = []
            else:
                m.append(line)

    if not len(smodels):
        raise RuntimeError("Could not extract any models from: {0}".format(dfile))

    extracted_models = []
    for i, m in enumerate(smodels):
        # TODO: Maybe change the name from quark to something a little more general
        fpath = os.path.join(directory, "quark_{0}.pdb".format(i))
        with open(fpath, 'w') as f:
            for l in m:
                # TODO: Reconstruct something sensible as from the coordinates on it's all quark-specific
                if l.startswith("ATOM"):
                    l = l[:54]+"  1.00  0.00              \n"
                f.write(l)
            extracted_models.append(fpath)
            LOGGER.debug("Wrote: {0}".format(fpath))
        
    return extracted_models

def tmp_file_name(delete=True, directory=None):
    """
    Return a filename for a temporary file

    Parameters
    ---------
    delete : bool
       Flag whether the temporary file should be deleted
    directory : str
       Path to a directory to write the files to.
    """
    directory = os.getcwd() if not directory else directory
    t = tempfile.NamedTemporaryFile(dir=directory, delete=delete)
    tmp1 = t.name
    t.close()
    return tmp1


# Header string
header ="""#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE - Ab Initio Modelling Molecular Replacement               #
#########################################################################

The authors of specific programs should be referenced where applicable:""" + \
"\n\n" + construct_references() + "\n\n"

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