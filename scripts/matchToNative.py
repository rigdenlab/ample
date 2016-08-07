#!/usr/bin/env ccp4-python

import logging
import os
import shutil
import sys

if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
sys.path.insert(0, os.path.join(os.environ['CCP4'], "share", "ample", "python"))
from ample.util import ample_util
from ample.util import csymmatch
from ample.util import mtz_util
from ample.util import pdb_edit
from ample.util import phenixer
from ample.util import shelxe

def run(nativePdb, nativeMtz=None, nativeMap=None, mrPdbs=None, outDir=None):
    
    logging.info("Using native PDB file: {0}".format(os.path.abspath(nativePdb)))
    logging.info("Using native MTZ file: {0}".format(os.path.abspath(nativeMtz)))
    logging.info("We will be wrapping the following MTZ files to the native: {0}".format(mrPdbs))
    
    phenix = False

    # Find out where we're running from
    if outDir is not None:
        if not os.path.isdir(outDir):
            raise RuntimeError, "Cannot find output directory: {0}".format(outDir)
        outDir = os.path.abspath(outDir)
    else:
        outDir = os.getcwd()
        
    if phenix:
        if nativeMap is None:
            nativeMap = generateMap(nativePdb, nativeMtz)
        if not os.path.isfile(nativeMap):
            raise RuntimeError, "Cannot find nativeMap: {0}".format(nativeMap)
    else:
        shelxeExe = ample_util.find_exe('shelxe')
    
    for mrPdb in mrPdbs:
        if phenix:
            logging.debug("Searching for origin shift using: {0} {1}".format(nativeMap, mrPdb))
            origin = phenixer.ccmtzOrigin(nativeMap, mrPdb)
            # offset.pdb is the mrPdb moved onto the new origin
            offsetPdb = "offset.pdb"
            logging.debug("Found origin: {0}\nOffset pdb is: {1}".format(origin, offsetPdb))
        else:
            originShift = shelxe.shelxe_origin(shelxeExe, nativePdb, nativeMtz, mrPdb)
            logging.debug("Found origin: {0}".format(originShift))
            offsetPdb = ample_util.filename_append(mrPdb, astr='offset', directory=os.getcwd())
            pdb_edit.translate(mrPdb, offsetPdb, originShift)
        
        # Run csymmatch to map offsetted pdb onto native
        csymmPdb = ample_util.filename_append(filename=mrPdb, astr="csymmatch", directory=outDir)
        logging.debug("Running csymmatch to wrap {0} onto native {1}".format(offsetPdb, nativePdb))
        csymmatch.Csymmatch().run(refPdb=nativePdb, inPdb=offsetPdb, outPdb=csymmPdb, originHand=False)
        
        logging.info("Matched PDB is: {0}".format(csymmPdb))
    
    return

def generateMap(nativePdb, nativeMtz):
    
    # Get the labels from the MTZ file
    logging.debug("Parsing MTZ file {0} to determine column labels".format(nativeMtz))
    F, SIGF, FREE = mtz_util.get_labels(nativeMtz)
    # Generate map from 
    logging.debug("Generating map from: {0} {1}".format(nativeMtz,
                                                 nativePdb))
    return phenixer.generateMap(nativeMtz,
                                nativePdb,
                                FP=F,
                                SIGFP=SIGF,
                                FREE=FREE)
    
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('native_pdb', type=str, help='Native PDB file')
    p.add_argument('native_mtz', type=str, help='Native MTZ file')
    p.add_argument('mr_pdbs', type=str, nargs='+', help='One or more PDB file from Molecular Replacement')
    args=p.parse_args()
    
    # Setup logging
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    
    run(nativePdb=args.native_pdb, nativeMtz=args.native_mtz, mrPdbs=args.mr_pdbs)

