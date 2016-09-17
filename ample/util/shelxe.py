#!/usr/bin/env ccp4-python
'''
Created on 2 Feb 2015

@author: jmht
'''
import logging
import os
import shutil
import sys

# our imports
from ample.util import ample_util
from ample.util import mtz_util

_logger = logging.getLogger(__name__)

if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
mrbumpd = os.path.join(os.environ['CCP4'],"share","mrbump","include","parsers")
sys.path.insert(0,mrbumpd)
import parse_shelxe

def shelxe_origin(shelxe_exe, native_pdb, native_mtz, mr_pdb):
    if not ample_util.is_exe(shelxe_exe): raise RuntimeError,"Cannot find shelxe executable: {0}".format(shelxe_exe)
    if not os.path.isfile(native_pdb): raise RuntimeError,"Cannot find native_pdb: {0}".format(native_pdb)
    if not os.path.isfile(native_mtz): raise RuntimeError,"Cannot find native_mtz: {0}".format(native_mtz)
    if not os.path.isfile(mr_pdb): raise RuntimeError,"Cannot find mr_pdb: {0}".format(mr_pdb)
    
    stem = "shelxe-input" # stem name for all shelxe files
    hkl_file = stem+".hkl"
    hkl_file = mtz_util.to_hkl(native_mtz, hkl_file=hkl_file)
    
    # Rename nativePdb and mrPdb
    shutil.copyfile(mr_pdb, stem+".pda")
    shutil.copyfile(native_pdb, stem+".ent")
    trace_cycles = 0
    frac_solvent = 0.5
    cmd = [shelxe_exe,'shelxe-input.pda','-a{0}'.format(trace_cycles),'-q', '-s{0}'.format(frac_solvent),'-o','-n','-t0','-m0','-x']
    logfile = os.path.abspath('shelxe.log')
    ret = ample_util.run_command(cmd=cmd, logfile=logfile, directory=None, dolog=True, stdin=None)
    if ret != 0:
        raise RuntimeError,"Error running shelxe - see log: {0}".format(logfile)
    else:
        for ext in ['.pda','.hkl','.ent','.pdo','.phs','.lst','_trace.ps']:
            try: os.unlink(stem+ext)
            except: pass
    
    sp = parse_shelxe.ShelxeLogParser(logfile)
    if not sp.originShift:
        raise RuntimeError,"SHELXE failed to find an origin. Please check the logfile: {0}".format(logfile)
    os.unlink(logfile)
    originShift=[ o*-1 for o in sp.originShift ]
    _logger.debug('shelxe_origin calculated origin: {0}'.format(originShift))
    return originShift

if __name__ == "__main__":

    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    import argparse
    parser = argparse.ArgumentParser(description='Determine origin using SHELXE', prefix_chars="-")
    
    parser.add_argument('--native_mtz',
                       help='Native MTZ', required=True)
    parser.add_argument('--native_pdb',
                       help='Native PDB', required=True)
    parser.add_argument('--mr_pdb',
                       help='Molecular Replacement MTZ', required=True)
    parser.add_argument('--executable',
                       help="Path to SHELXE executable")

    args = parser.parse_args()
    
    executable = None
    if args.executable:
        executable = args.executable
    else:
        executable = os.path.join(os.environ['CCP4'],"bin","shelxe" + ample_util.EXE_EXT)
    
    # Get full paths to all files
    native_mtz = os.path.abspath(args.native_mtz)
    if not os.path.isfile(native_mtz):
        raise RuntimeError, "Cannot find input file: {0}".format(native_mtz)
    native_pdb = os.path.abspath(args.native_pdb)
    if not os.path.isfile(native_pdb):
        raise RuntimeError, "Cannot find input file: {0}".format(native_pdb)
    mr_pdb = os.path.abspath(args.mr_pdb)
    if not os.path.isfile(mr_pdb):
        raise RuntimeError, "Cannot find input file: {0}".format(mr_pdb)
    
    origin_shift = shelxe_origin(executable, native_pdb, native_mtz, mr_pdb)
    print "Origin shift is: {0}".format(origin_shift)
        
