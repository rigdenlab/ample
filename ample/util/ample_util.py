"""Various miscellaneous functions"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Jan 2016"
__version__ = "1.0"

import cPickle
import glob
import logging
import os
import subprocess
import sys
import tarfile
import tempfile
import warnings
import zipfile

import exit_util
import pdb_edit

from ample.constants import SHARE_DIR

CCP4_VERSION=None
SCRIPT_EXT = '.bat' if sys.platform.startswith('win') else '.sh'
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'

class FileNotFoundError(Exception): pass

# ample_util is used before anything else so there is no logger available
# and we need to a Null handler
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def ccp4_version():
    """Get CCP4 version as a tuple

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


def construct_references(optd):
    """Construct the reference list

    Description
    -----------
    This is somewhat a very basic BibTex parser. It is under no circumstances foolproof but rather
    the bare minimum for what is required in AMPLE.

    If a more sophisticated parser is required, refer to one of the python packages available online.

    Parameters
    ----------
    optd : dict
       A dictionary containing AMPLE options

    Returns
    -------
    str
        A string containing the references

    """
    # ========================================
    # Get the filename and check we can use it
    # ========================================
    ref_fname = os.path.join(SHARE_DIR, "include", "ample.bib")
    if not is_file(ref_fname):
        msg = "Cannot find BibTex file containing references. " \
              "Please determine them yourself and cite AMPLE."
        return msg

    articles = []
    article = {}
    entry = False

    with open(ref_fname, "r") as fhin:
        for line in fhin.readlines():

            line = line.strip()

            if not line:                                # Make sure line is not empty
                continue

            elif line.startswith("@"):                  # Beginning of all BibTex entry blocks
                entry = True                            # Notify that we have an article block
                unique_id = line.replace("@article{", "").replace(",", "")
                article = {'unique_id': unique_id}      # Reset the article dictionary

            elif line == "}":                           # End of all BibTex entry blocks
                entry = False                           # Notify that we article block is over
                articles.append(article)                # Save the article dictionary

            elif entry:                                 # BibTex entry block
                # Some dirty line handling.
                # Not very bulletproof but should do for now
                line = line.replace("{", "").replace("}", "")
                key, value = [l.strip() for l in line.split("=")]
                value = value.rstrip(",").replace("\"", "")
                # Save the data to the article entry
                article[key] = value

    # =====================================================
    # Get the used applications to not print all references
    # =====================================================
    def used_apps(optd):
        """determine which applications were used and return their labels"""

        # For now print all AMPLE papers
        labels = ['AMPLE', 'AMPLE_MODELLING', 'AMPLE_COILS', 'AMPLE_CONTACTS', 'CCP4']

        # This conditional attempts to recognise which programs have run and therefore
        # prints only the relevant references. It is under no circumstance perfect.
        if optd['nmr_model_in']:
            labels += ['AMPLE_NMR']

        # Flags related to cluster-and-truncate approach
        elif not optd['import_ensembles']:
            labels += ['THESEUS']
            if optd['use_scwrl']:
                labels += ['SCWRL4']
            elif optd['cluster_method'] in ['spicker', 'spicker_qscore', 'spicker_tmscore']:
                labels = ['FPC']
            elif optd['cluster_method'] in ['fast_protein_cluster']:
                labels += ['SPICKER']

        # Flags related to Molecular Replacement
        elif optd['do_mr']:
            labels += ['MrBUMP']
            if optd['use_arpwarp']:
                labels += ['Arp_Warp']
            elif optd['use_buccaneer']:
                labels += ['Buccaneer']
            elif optd['use_shelxe']:
                labels += ['SHELXE']
            elif 'molrep' in optd['mrbump_programs']:
                labels += ['Molrep', 'REFMAC']
            elif 'phaser' in optd['mrbump_programs']:
                labels += ['Phaser', 'REFMAC']

        return labels

    applications = used_apps(optd)

    to_remove = [i for i, art in enumerate(articles) if art['label'] not in applications]
    # reverse order so that we don't throw off the subsequent indexes
    for index in sorted(to_remove, reverse=True):
        del articles[index]

    # =========================================================================
    # Somewhat a template of how we want to write each article in BibTex format
    # =========================================================================
    template_bib = "@article{{{unique_id},{sep}author = {{{author}}},{sep}doi = {{{doi}}},{sep}" \
                   "journal = {{{journal}}},{sep}number = {{{number}}},{sep}pages = {{{pages}}},{sep}" \
                   "title = {{{{{title}}}}},{sep}volume = {{{volume}}},{sep}year = {{{year}}},{sep}}}{sep}"
    references_bib = [template_bib.format(sep=os.linesep, **article) for article in articles]

    ref_fname = os.path.join(optd['work_dir'], optd['name']+".bib")
    with open(ref_fname, "w") as fhout:
        fhout.write(os.linesep.join(references_bib))

    # ==========================================================
    # Somewhat a template of how we want to display each article
    # ==========================================================
    for article in articles:
        # Shorten the author list for the displayed message
        article['author'] = article['author'].split(" and ")[0].split(",")[0] + " et al."
        # Display page separator as single dash
        article['pages'] = article['pages'].replace("--", "-")
    template_msg = "* {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
    references_msg = [template_msg.format(**article) for article in articles]

    return (os.linesep*2).join(references_msg)


def extract_models(amoptd, sequence=None, single=True, allsame=True):
    """Extract some models
    
    Description
    -----------
    Check a directory of pdbs or extract pdb files from a given tar/zip file or directory of pdbs
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
            msg = "Cannot find models file: {0}".format(filename)
            exit_util.exit_error(msg)
            
        # we need a models_dir to extract into
        assert models_dir, "extract_models() needs a models_dir path!"
        if not os.path.isdir(models_dir):
            os.mkdir(models_dir)
        models_dir = models_dir
        
        # See what sort of file this is:
        f, suffix = os.path.splitext(filename)
        if suffix in ['.gz', '.bz']:
            f, s2 = os.path.splitext(f)
            if s2 == '.tar':
                suffix = s2+suffix
        
        tsuffixes = ['.tar.gz', '.tgz', '.tar.bz', '.tbz']
        suffixes = tsuffixes + ['.zip']
        if suffix not in suffixes:
            msg = "Do not know how to extract files from file: {0}\n " \
                  "Acceptable file types are: {1}".format(filename, suffixes)
            exit_util.exit_error(msg)
        if suffix in tsuffixes:
            files = extract_tar(filename, models_dir)
        else:
            files = extract_zip(filename, models_dir)
        
        # Assume anything with one member is quark decoys
        if len(files) == 1:
            quark_filename = 'alldecoy.pdb'
            f = os.path.basename(files[0])
            if not f == quark_filename:
                msg = "Only found one member ({0}) in file: {1} " \
                      "and the name was not {2}\n".format(f, filename, quark_filename)
                msg += "If this file contains valid QUARK decoys, please email: ccp4@stfc.ac.uk"
                exit_util.exit_error(msg)
            # Now extract the quark pdb files from the monolithic file
            split_quark(files[0], models_dir)
            # We delete the quark_name file as otherwise we'll try and model it
            os.unlink(files[0])
            # If we've got quark models we don't want to modify the side chains as we only have polyalanine so we
            # set this here - horribly untidy as we should have one place to decide on side chains
            logger.info('Found QUARK models in file: {0}'.format(filename))
            amoptd['quark_models'] = True
    
    if not pdb_edit.check_pdb_directory(models_dir, sequence=sequence, single=single, allsame=allsame):
        msg = "Problem importing pdb files - please check the log for more information"
        exit_util.exit_error(msg)
    
    amoptd['models_dir'] = models_dir
    return glob.glob(os.path.join(models_dir, "*.pdb"))


def extract_tar(filename, directory, suffixes=['.pdb']):
    # Extracting tarfile
    logger.info('Extracting files from tarfile: {0}'.format(filename) )
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
    logger.info('Extracting files from zipfile: {0}'.format(filename) )
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
    
    Parameters
    ----------
    executable : str
       The name of the program or the path to an existing executable
    dirs : list, tuple, optional
       Additional directories to search for the location
    """
    logger.debug('Looking for executable: {0}'.format(executable) )
    
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
        logger.debug('Checking paths: {0}'.format(paths))
        
        for path in paths:
            exe_file = os.path.abspath(os.path.join(path, executable))   
            if is_exe(exe_file):
                logger.debug( 'Found executable {0} in directory {1}'.format(executable,path) )
                found=True
                break
    
    if not found:
        raise FileNotFoundError("Cannot find executable: {0}".format(executable))
    
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


def ideal_helices(optd):
    """Get some ideal helices

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
    nresidues = optd['fasta_length']
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
        
    optd['ensembles'] = pdbs
    optd['ensemble_options'] = ensemble_options
    optd['ensembles_data'] = ensembles_data
    return


def is_exe(fpath):
    """Check if an executable exists

    Parameters
    ----------
    fpath : str
       The path to the executable
    
    Returns
    -------
    bool

    """
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)


def is_file(fpath):
    """Check if a file exists

    Parameters
    ----------
    fpath : str
       The path to the file

    Returns
    ------
    bool

    """
    return fpath and os.path.isfile(fpath) and \
           os.access(fpath, os.R_OK) and os.stat(fpath).st_size > 0


def make_workdir(work_dir, ccp4_jobid=None, rootname='AMPLE_'):
    """Make a work directory rooted at work_dir and return its path

    Parameters
    ----------
    work_dir : str
       The path to a working directory
    ccp4_jobid : int, optional
       CCP4-assigned job identifier
    rootname : str, optional
        Base name of the AMPLE directory [default: \'AMPLE_\']

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
    """Execute a command and return the exit code.

    Parameters
    ----------
    cmd : list
       Command to run as a list
    stdin : str, optional
       Stdin for the command
    logfile : str, optional
       The path to the logfile
    directory : str, optional
       The directory to run the job in (cwd assumed)
    dolog : bool, optional
       Whether to output info to the system log [default: False]

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
        logger.debug("In directory {0}".format(directory))
        logger.debug("Running command: {0}".format(" ".join(cmd)))
        if kwargs:
            logger.debug("kwargs are: {0}".format(kwargs))

    file_handle = False
    if logfile:
        if type(logfile) == file:
            file_handle = True
            logf = logfile
            logfile = os.path.abspath(logf.name)
        else:
            logfile = os.path.abspath(logfile)
            logf = open(logfile, "w")
        if dolog: logger.debug("Logfile is: {0}".format(logfile))
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
        if dolog: logger.debug("stdin for cmd was: {0}".format(stdinstr))

    p.wait()
    if not file_handle:
        logf.close()

    return p.returncode


def read_amoptd(amoptd_fname):
    """Read a PICKLE-formatted AMPLE options file

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
        raise RuntimeError("Cannot access AMPLE options file: {0}\n".format(amoptd_fname))

    with open(amoptd_fname, 'r') as f:
        amoptd = cPickle.load(f)
        logger.info("Loaded state from file: {0}\n".format(amoptd['results_path']))
    return amoptd


def saveAmoptd(*args):
    """Save AMPLE options to a PICKLE-formatted file

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
    """Save AMPLE options to a PICKLE-formatted file

    Parameters
    ----------
    amoptd : dict
       AMPLE options from saved state

    """
    # Save results
    with open(amoptd['results_path'], 'w') as f:
        cPickle.dump(amoptd, f)
        logger.info("Saved state as file: {0}\n".format(amoptd['results_path']))
    return


def split_quark(*args):
    """Split a single PDB with multiple models in individual PDB files

    See Also
    --------
    split_models

    """
    return split_models(*args)


def split_models(dfile, directory):
    """Split a single PDB with multiple models in individual PDB files

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
    logger.info("Extracting decoys from: {0} into {1}".format(dfile, directory))
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
            logger.debug("Wrote: {0}".format(fpath))
        
    return extracted_models


def tmpFileName():
    """Return a filename for a temporary file

    See Also
    --------
    tmp_file_name

    Warnings
    --------
    This function was deprecated and will be removed in future releases. Please use ``tmp_file_name()`` instead.

    """
    msg = "This function was deprecated and will be removed in future release"
    warnings.warn(msg, DeprecationWarning, stacklevel=2)
    return tmp_file_name()


def tmp_file_name(delete=True, directory=None, suffix=""):
    """Return a filename for a temporary file

    Parameters
    ----------
    delete : bool, optional
       Flag whether the temporary file should be deleted [default: True]
    directory : str, optional
       Path to a directory to write the files to.
    suffix : str, optional
       A suffix to the temporary filename

    """
    directory = os.getcwd() if not directory else directory
    t = tempfile.NamedTemporaryFile(dir=directory, delete=delete, suffix=suffix)
    tmp1 = t.name
    t.close()
    return tmp1

# ======================================================================
# Some default string messages that we need during the program to inform
# the user of certain information
# ======================================================================

header = """
#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE - Ab Initio Modelling Molecular Replacement               #
#########################################################################

"""

# ======================================================================
# ======================================================================

reference = """
#########################################################################

The authors of specific programs should be referenced where applicable:

{refs}

"""

# ======================================================================
# ======================================================================

survey_url = "http://goo.gl/forms/7xP9M4P81O"
footer = """

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
# ======================================================================
# ======================================================================
