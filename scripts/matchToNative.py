#!/usr/bin/env ccp4-python

import logging
import os
import shutil
import sys

if not "CCP4" in os.environ.keys():
    raise RuntimeError('CCP4 not found')
sys.path.insert(0, os.path.join(os.environ['CCP4'], "share",
                                "ample", "python"))

from ample.util import ample_util
from ample.util import csymmatch
from ample.util import mtz_util
from ample.util import pdb_edit
from ample.util import phenixer
from ample.util import shelxe


def run(nativePdb, nativeMtz=None, nativeMap=None, mrPdbs=None, outDir=None):

    logging.info("Using native PDB file: %s", os.path.abspath(nativePdb))
    logging.info("Using native MTZ file: %s", os.path.abspath(nativeMtz))
    logging.info("Wrapping the following MTZ files to the native: %s",
                 " ".join(mrPdbs))

    phenix = False

    if outDir is not None:
        if not os.path.isdir(outDir):
            msg = "Cannot find output directory: {}".format(outDir)
            raise RuntimeError(msg)
        outDir = os.path.abspath(outDir)
    else:
        outDir = os.getcwd()

    if phenix:
        if nativeMap is None:
            nativeMap = generateMap(nativePdb, nativeMtz)
        if not os.path.isfile(nativeMap):
            msg = "Cannot find nativeMap: {}".format(nativeMap)
            raise RuntimeError(msg)
    else:
        shelxeExe = ample_util.find_exe('shelxe')

    removables = []
    for mrPdb in mrPdbs:
        if phenix:
            logging.debug("Searching for origin shift using: %s %s",
                          nativeMap, mrPdb)
            origin = phenixer.ccmtzOrigin(nativeMap, mrPdb)
            # offset.pdb is the mrPdb moved onto the new origin
            offsetPdb = "offset.pdb"
            logging.debug("Found origin: %s\nOffset pdb is: %s",
                          origin, offsetPdb)
        else:

            mrinfo = shelxe.MRinfo(shelxeExe, nativePdb, nativeMtz)
            mrinfo.analyse(mrPdb)
            originShift = mrinfo.originShift
            logging.debug("Found origin: {0}".format(originShift))
            offsetPdb = ample_util.filename_append(mrPdb, astr='offset',
                                                   directory=os.getcwd())
            pdb_edit.translate(mrPdb, offsetPdb, originShift)

        csymmPdb = ample_util.filename_append(filename=mrPdb, astr="csymmatch",
                                              directory=outDir)
        logging.debug("Running csymmatch to wrap %s onto native %s",
                      offsetPdb, nativePdb)
        csymmatch.Csymmatch().run(refPdb=nativePdb, inPdb=offsetPdb,
                                  outPdb=csymmPdb, originHand=False)
        removables += [offsetPdb]
        logging.info("Matched PDB is: %s", csymmPdb)

    map(os.remove, removables + ["shelxe-input.hkl", "shelxe-input.ent"])


def generateMap(nativePdb, nativeMtz):
    logging.debug("Parsing MTZ file %s to determine column labels", nativeMtz)
    F, SIGF, FREE = mtz_util.get_labels(nativeMtz)
    logging.debug("Generating map from: %s %s", nativeMtz, nativePdb)
    return phenixer.generateMap(nativeMtz, nativePdb, FP=F, SIGFP=SIGF,
                                FREE=FREE)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-odir', type=str, help="output directory")
    p.add_argument('native_pdb', type=str, help='Native PDB file')
    p.add_argument('native_mtz', type=str, help='Native MTZ file')
    p.add_argument('mr_pdbs', type=str, nargs='+',
                   help='One or more PDB file from Molecular Replacement')
    args = p.parse_args()

    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    run(nativePdb=args.native_pdb, nativeMtz=args.native_mtz,
        mrPdbs=args.mr_pdbs, outDir=args.odir)
