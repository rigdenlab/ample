#!/usr/bin/env ccp4-python
'''
Created on 2 Feb 2015

@author: jmht
'''
import os
import shutil
import sys
import uuid

from ample.util import ample_util
from ample.util import mtz_util


try:
    from mrbump.parsers import parse_shelxe
except ImportError:
    mrbumpd = os.path.join(os.environ['CCP4'], "share", "mrbump", "include", "parsers")
    sys.path.insert(0,mrbumpd)
    import parse_shelxe


class MRinfo(object):
    """An object to analyse Molecular Replacement solutions

    Attributes
    ----------
    work_dir : str
      Path to the working directory
    shelxe_exe : str
      Path to the SHELXE executable
    stem : str
      The name for all the SHELXE files
    originShift : list
      The origin shift of the MR pdb to native as a list of three floats
    MPE : float
      The Mean Phase Error of the MR pdb to the native pdb
    wMPE : float
      The weighted Mean Phase Error of the MR pdb to the native pdb

    """
    def __init__(self, shelxe_exe, native_pdb, native_mtz, work_dir=None):
        """Intialise from native pdb and mtz so that analyse only requires a MR pdb

        Parameters
        ----------
        shelxe_exe : str
          Path to the SHELXE executable
        native_pdb : str
          Path to the native PDB file
        native_mtz : str
          Path to the native MTZ file

        """
        if work_dir is None: work_dir = os.getcwd()
        self.work_dir = work_dir
        self.shelxe_exe = shelxe_exe
        self.stem = 'shelxe-input-{}'.format(str(uuid.uuid1()))

        self.MPE = None
        self.wMPE = None
        self.originShift = None

        self.mk_native_files(native_pdb, native_mtz)
        return

    def mk_native_files(self, native_pdb, native_mtz):
        """Create the files required by SHELXE from the native structure

        Parameters
        ----------
        native_pdb : str
          Path to the native PDB file
        native_mtz : str
          Path to the native MTZ file

        """
        mtz_util.to_hkl(native_mtz, hkl_file=os.path.join(self.work_dir, self.stem + ".hkl"))
        shutil.copyfile(native_pdb, os.path.join(self.work_dir, self.stem + ".ent"))

    def analyse(self, mr_pdb, cleanup=True):
        """Use SHELXE to analyse an MR pdb file to determine the origin shift and phase error

        This function sets the ``MPE``, ``wMPE`` and ``originShift`` attributes.

        Parameters
        ----------
        mr_pdb : str
          Path to the Molecular Replacement PDB file

        """

        os.chdir(self.work_dir)
        input_pdb = self.stem + ".pda"
        shutil.copyfile(mr_pdb, os.path.join(self.work_dir, input_pdb))

        cmd = [self.shelxe_exe, input_pdb, '-a0', '-q', '-s0.5', '-o', '-n', '-t0', '-m0', '-x']
        logfile = os.path.abspath('shelxe_{}.log'.format( str(uuid.uuid1())))
        ret = ample_util.run_command(cmd=cmd, logfile=logfile, directory=None, dolog=False, stdin=None)
        if ret != 0: 
            raise RuntimeError("Error running shelxe - see log: {0}".format(logfile))

        sp = parse_shelxe.ShelxeLogParser(logfile)
        # Only added in later version of MRBUMP shelxe parser
        if hasattr(sp, 'MPE'):
            self.MPE = sp.MPE
        self.wMPE = sp.wMPE
        if isinstance(sp.originShift, list):
            self.originShift = [ o*-1 for o in sp.originShift ]

        if cleanup:
            for ext in ['.hkl', '.ent', '.pda','.pdo','.phs','.lst','_trace.ps']:
                try:
                    os.unlink(self.stem + ext)
                except:
                    pass
            os.unlink(logfile)

if __name__ == "__main__":

    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    import argparse
    parser = argparse.ArgumentParser(description='Determine origin using SHELXE', prefix_chars="-")
    parser.add_argument('--native_mtz', help='Native MTZ', required=True)
    parser.add_argument('--native_pdb', help='Native PDB', required=True)
    parser.add_argument('--mr_pdb', help='Molecular Replacement MTZ', required=True)
    parser.add_argument('--executable', help="Path to SHELXE executable")
    args = parser.parse_args()

    executable = None
    if args.executable:
        executable = args.executable
    else:
        executable = os.path.join(os.environ['CCP4'],"bin","shelxe" + ample_util.EXE_EXT)

    # Get full paths to all files
    native_mtz = os.path.abspath(args.native_mtz)
    if not os.path.isfile(native_mtz):
        raise RuntimeError("Cannot find input file: {0}".format(native_mtz))
    native_pdb = os.path.abspath(args.native_pdb)
    if not os.path.isfile(native_pdb):
        raise RuntimeError("Cannot find input file: {0}".format(native_pdb))
    mr_pdb = os.path.abspath(args.mr_pdb)
    if not os.path.isfile(mr_pdb):
        raise RuntimeError("Cannot find input file: {0}".format(mr_pdb))

    mrinfo = MRinfo(executable, native_pdb, native_mtz)
    mrinfo.analyse(mr_pdb)
    os.unlink('shelxe-input.hkl')
    os.unlink('shelxe-input.ent')
    print("Origin shift is: {0}".format(mrinfo.originShift))
