"""Various miscellaneous functions"""
__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Jan 2016"
__version__ = "1.0"

from contextlib import contextmanager
import pickle
import logging
import os
import subprocess
import sys
import tarfile
import tempfile
import warnings
import zipfile

from ample.util import ccp4, exit_util

from ample.constants import SHARE_DIR, AMPLEDIR, I2DIR

CCP4 = ccp4.CCP4()
SCRIPT_EXT = '.bat' if sys.platform.startswith('win') else '.sh'
EXE_EXT = '.exe' if sys.platform.startswith('win') else ''
SCRIPT_HEADER = '' if sys.platform.startswith('win') else '#!/bin/bash'


class FileNotFoundError(Exception):
    pass


# ample_util is used before anything else so there is no logger available
# and we need to a Null handler
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

three2one = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK': 'X',
}
# http://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
# aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )
one2three = dict((v, k) for (k, v) in three2one.items())


def amoptd_fix_path(optd, newroot):
    """Update all the paths in an AMPLE results dictionary to be rooted at newroot

    Parameters
    ----------
    optd: dict
      AMPLE results dictionary
    newroot: str
       Path to the AMPLE root directory (topdir containing MRBUMP dir etc.)
    """
    # oldroot = os.sep.join(optd['work_dir'].split(os.sep)[:-1])
    oldroot = optd['work_dir']
    for k in ['benchmark_dir', 'native_pdb', 'native_pdb_std', 'fasta']:
        if k in optd and isinstance(optd[k], str):
            optd[k] = optd[k].replace(oldroot, newroot)

    MRBUMP_FILE_KEYS = [
        'PHASER_logfile',
        'PHASER_pdbout',
        'PHASER_mtzout',
        'REFMAC_logfile',
        'REFMAC_pdbout',
        'REFMAC_mtzout',
        'BUCC_logfile',
        'BUCC_pdbout',
        'BUCC_mtzout',
        'ARP_logfile',
        'ARP_pdbout',
        'ARP_mtzout',
        'SHELXE_logfile',
        'SHELXE_pdbout',
        'SHELXE_mtzout',
        'SXRBUCC_logfile',
        'SXRBUCC_pdbout',
        'SXRBUCC_mtzout',
        'SXRARP_logfile',
        'SXRARP_pdbout',
        'SXRARP_mtzout',
    ]
    if 'mrbump_results' in optd:
        for r in optd['mrbump_results']:
            for k in MRBUMP_FILE_KEYS:
                if k in r and isinstance(r[k], str):
                    warnings.warn("FIX MRBUMP BUG buccaneer refine.pdb vs refined.pdb")
                    new = r[k].replace(oldroot, newroot)
                    # logger.info('Changing amopt entry %s from: %s to: %s', k, old, new)
                    r[k] = new
    return optd


@contextmanager
def disable_logging(logger, max_loglevel=logging.CRITICAL):
    """A context manager to disable logging within a block.
    
    Parameters
    ----------
    logger : logging.Logger
       logger instance
    highet_level: int
        loglevel that will be set for the duration of the context
    """
    previous_level = None
    if logger.getEffectiveLevel() < max_loglevel:
        previous_level = logger.level
        logger.setLevel(max_loglevel)
    try:
        yield
    finally:
        if previous_level is not None:
            # changed loglevel so reset it
            logger.setLevel(previous_level)


def extract_tar(archive, directory=None, filenames=None, suffixes=None):
    """Extract one or more files from a tar file into a specified directory

    Parameters
    ----------
    archive : str
       tar archive to extract from
    directory : str
       directory to extract files into
    filenames : list
       a list of files to extract from the archive
    suffixes : list
       only extract a file if the suffix is in the list


    Returns
    -------
    list
        A list of the extracted files
    """

    def extract_me(member, filenames, suffixes):
        """If filenames or suffixes is given, only extract files with those filenames or suffixes,
        otherwise extract all files."""
        if filenames:
            if os.path.basename(member.name) in filenames:
                return True
        else:
            if suffixes:
                if os.path.splitext(member.name)[1] in suffixes:
                    return True
            else:
                return True
        return False

    if not directory:
        directory = os.getcwd()
    if not os.path.isdir(directory):
        os.mkdir(directory)

    logger.info('Extracting files from tarfile: %s into directory: %s', archive, directory)
    files = []
    with tarfile.open(archive, 'r:*') as tf:
        members = tf.getmembers()
        if members:
            for member in members:
                if extract_me(member, filenames, suffixes):
                    member.name = os.path.basename(member.name)  # Hack to remove any paths
                    tf.extract(member, path=directory)
                    files.append(os.path.abspath(os.path.join(directory, member.name)))
        else:
            logger.critical('Empty archive: %s', archive)
    return files


def extract_zip(filename, directory, suffixes=None):
    # zip file extraction
    logger.info('Extracting files from zipfile: %s', filename)
    if not zipfile.is_zipfile(filename):
        msg = 'File is not a valid zip archive: {0}'.format(filename)
        exit_util.exit_error(msg)
    zipf = zipfile.ZipFile(filename)
    zif = zipf.infolist()
    if not zif:
        msg = 'Empty zip file: {0}'.format(filename)
        exit_util.exit_error(msg)
    files = []
    for f in zif:
        if os.path.splitext(f.filename)[1] in suffixes:
            # Hack to rewrite name
            f.filename = os.path.basename(f.filename)
            zipf.extract(f, path=directory)
            files.append(os.path.join(directory, f.filename))
    if not files:
        msg = 'Could not find any files with suffixes {0} in zipfile: {1}'.format(suffixes, filename)
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
    logger.debug('Looking for executable: %s', executable)

    exe_file = None
    found = False
    if is_exe(executable):
        exe_file = os.path.abspath(executable)
        found = True
    else:
        # If the user has given a path we just take the name part
        _, fname = os.path.split(executable)
        if fname:
            executable = fname

        # By default we search in the system PATH and add any additional user given paths here
        paths = os.environ["PATH"].split(os.pathsep)
        if dirs:
            paths += dirs
        logger.debug('Checking paths: %s', paths)

        for path in paths:
            exe_file = os.path.abspath(os.path.join(path, executable))
            if is_exe(exe_file):
                logger.debug('Found executable %s in directory %s', executable, path)
                found = True
                break

    if not found:
        raise FileNotFoundError("Cannot find executable: {0}".format(executable))

    logger.debug('find_exe found executable: %s', exe_file)
    return exe_file


def filename_append(filename=None, astr=None, directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split(filename)
    name, suffix = os.path.splitext(fname)
    name = name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join(directory, name)


def ideal_helices(optd, ensembles=False):
    """Get some ideal helices

    Parameters
    ----------
    nresidues : int
       Number of residues to be used

    ensembles : bool
       If True then helical ensembles will be used instead of polyalanine ideal helices

    Returns
    -------
    pdbs : list
    ensemble_options : dict
    ensembles_data : list
    """
    nresidues = optd['fasta_length']
    include_dir = os.path.join(SHARE_DIR, 'include')
    if not ensembles:
        names = ['polyala_5', 'polyala_10', 'polyala_15', 'polyala_20', 'polyala_25', 'polyala_30', 'polyala_35',
                 'polyala_40']
    else:
        names = ['ensemble_20_gradientbfact_homogenous', 'ensemble_20_gradientbfact_heterogenous',
                 'ensemble_20_nativebfact_heterogenous', 'ensemble_20_nativebfact_homogenous',
                 'ensemble_20_rmsd_perpositionbfact_heterogenous', 'ensemble_20_rmsd_permodelbfact_homogenous',
                 'ensemble_20_rmsd_permodelbfact_heterogenous', 'ensemble_20_rmsd_perpositionbfact_homogenous',
                 'ensemble_15_gradientbfact_heterogenous', 'ensemble_15_gradientbfact_homogenous',
                 'ensemble_15_nativebfact_heterogenous', 'ensemble_15_nativebfact_homogenous',
                 'ensemble_15_rmsd_perpositionbfact_heterogenous',
                 'ensemble_15_rmsd_permodelbfact_heterogenous', 'ensemble_15_rmsd_perpositionbfact_homogenous',
                 'ensemble_15_rmsd_permodelbfact_homogenous', 'ensemble_25_gradientbfact_homogenous',
                 'ensemble_25_gradientbfact_heterogenous', 'ensemble_25_nativebfact_heterogenous',
                 'ensemble_25_nativebfact_homogenous', 'ensemble_25_rmsd_perpositionbfact_homogenous',
                 'ensemble_25_rmsd_permodelbfact_homogenous', 'ensemble_25_rmsd_perpositionbfact_heterogenous',
                 'ensemble_25_rmsd_permodelbfact_heterogenous', 'ensemble_30_gradientbfact_homogenous',
                 'ensemble_30_gradientbfact_heterogenous', 'ensemble_30_nativebfact_heterogenous',
                 'ensemble_30_nativebfact_homogenous', 'ensemble_30_rmsd_perpositionbfact_homogenous',
                 'ensemble_30_rmsd_perpositionbfact_heterogenous',
                 'ensemble_30_rmsd_permodelbfact_heterogenous', 'ensemble_30_rmsd_permodelbfact_homogenous',
                 'ensemble_35_gradientbfact_homogenous', 'ensemble_35_gradientbfact_heterogenous',
                 'ensemble_35_nativebfact_heterogenous', 'ensemble_35_nativebfact_homogenous',
                 'ensemble_35_rmsd_perpositionbfact_heterogenous',
                 'ensemble_35_rmsd_perpositionbfact_homogenous', 'ensemble_35_rmsd_permodelbfact_homogenous',
                 'ensemble_35_rmsd_permodelbfact_heterogenous', 'ensemble_40_gradientbfact_homogenous',
                 'ensemble_40_gradientbfact_heterogenous', 'ensemble_40_nativebfact_heterogenous',
                 'ensemble_40_nativebfact_homogenous', 'ensemble_40_rmsd_perpositionbfact_heterogenous',
                 'ensemble_40_rmsd_permodelbfact_heterogenous', 'ensemble_40_rmsd_perpositionbfact_homogenous',
                 'ensemble_40_rmsd_permodelbfact_homogenous', 'ensemble_10_gradientbfact_heterogenous',
                 'ensemble_10_gradientbfact_homogenous', 'ensemble_10_nativebfact_heterogenous',
                 'ensemble_10_nativebfact_homogenous', 'ensemble_10_rmsd_permodelbfact_homogenous',
                 'ensemble_10_rmsd_perpositionbfact_homogenous', 'ensemble_10_rmsd_permodelbfact_heterogenous',
                 'ensemble_10_rmsd_perpositionbfact_heterogenous', 'ensemble_5_gradientbfact_heterogenous',
                 'ensemble_5_gradientbfact_homogenous', 'ensemble_5_nativebfact_heterogenous',
                 'ensemble_5_nativebfact_homogenous', 'ensemble_5_rmsd_perpositionbfact_homogenous',
                 'ensemble_5_rmsd_permodelbfact_heterogenous', 'ensemble_5_rmsd_permodelbfact_homogenous',
                 'ensemble_5_rmsd_perpositionbfact_heterogenous']

    ensemble_options = {}
    ensembles_data = []
    pdbs = []
    for name in names:
        nres = int(name.split('_')[1])
        ncopies = nresidues / nres
        if ncopies < 1:
            ncopies = 1
        ensemble_options[name] = {'ncopies': ncopies}
        pdb = os.path.join(include_dir, "{0}.pdb".format(name))
        # Needed for pyrvapi results
        ensembles_data.append({'name': name, 'ensemble_pdb': pdb, 'num_residues': nres})
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
    return fpath and os.path.isfile(fpath) and os.access(fpath, os.R_OK) and os.stat(fpath).st_size > 0


def make_workdir(run_dir, ccp4i2=False, MAX_WORKDIRS=100):
    """Make a work directory rooted at run_dir and return its path

    Parameters
    ----------
    run_dir : str
       The path to a run directory where the job was started
    ccp4i2 : bool, optional
        Indicate if we are running under CCP4I2

    Returns
    -------
    work_dir : str
       The path to the working directory

    """
    if ccp4i2:
        work_dir = os.path.join(run_dir, I2DIR)
    else:
        run_inc = 0
        while True:
            work_dir = os.path.join(run_dir, AMPLEDIR + str(run_inc))
            if not os.path.exists(work_dir):
                break
            run_inc += 1
            if run_inc > MAX_WORKDIRS:
                raise RuntimeError("Too many work directories! {0}".format(work_dir))
    if os.path.exists(work_dir):
        raise RuntimeError(
            "There is an existing AMPLE work directory: {0}\n" "Please delete/move it aside.".format(work_dir)
        )
    os.mkdir(work_dir)
    return work_dir


def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False, **kwargs):
    """Execute a command and return the exit code.

    Parameters
    ----------
    cmd : list
       Command to run as a list
    stdin : str or filehandle, optional
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
        logger.debug("In directory %s", directory)
        logger.debug("Running command: %s", " ".join(cmd))
        if kwargs:
            logger.debug("kwargs are: %s", str(kwargs))

    file_handle = False
    if logfile:
        try:
            logfile = os.path.abspath(logfile)
        except TypeError:
            file_handle = True
        if file_handle:
            logf = logfile
            logfile = os.path.abspath(logf.name)
        else:
            logf = open(logfile, "w")
            logfile = os.path.abspath(logfile)
        if dolog:
            logger.debug("Logfile is: %s", logfile)
    else:
        logf = tempfile.NamedTemporaryFile(dir=directory, delete=False, suffix='')
    if stdin is not None:
        stdinstr = stdin
        stdin = subprocess.PIPE
    # Windows needs some special treatment
    if os.name == "nt":
        kwargs.update({'bufsize': 0, 'shell': "False"})
    p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs)

    if stdin is not None:
        p.stdin.write(stdinstr.encode())
        p.stdin.close()
        if dolog:
            logger.debug("stdin for cmd was: %s", stdinstr)
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
        amoptd = pickle.load(f)
        logger.info("Loaded state from file: %s\n", amoptd['results_path'])
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
        pickle.dump(amoptd, f)
        logger.info("Saved state as file: %s\n", amoptd['results_path'])


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
