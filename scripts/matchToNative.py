#!/usr/bin/env ccp4-python

import os
import shutil
import sys

root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
sys.path.insert( 0, os.path.join( root, "python" ) )
sys.path.insert( 0, os.path.join( root, "scripts" ) )


import ample_util
import csymmatch
import mtz_util
import pdb_edit
import phenixer
import shelxe

# # We use MRBUMP's MTZ_parse
# sys.path.append(os.path.join(os.environ["CCP4"], "share", "mrbump", "include", "file_info")) # For MTZ_parse
# import MTZ_parse

def run( nativePdb, nativeMtz=None, nativeMap=None, mrPdbs=None, outDir=None ):
    
    phenix=False

    # Find out where we're running from
    if outDir is not None:
        if not os.path.isdir(outDir):
            raise RuntimeError,"Cannot find output directory: {0}".format(outDir)
        outDir=os.path.abspath(outDir)
    else:
        outDir=os.getcwd()
        
    if phenix:
        if nativeMap is None:
            nativeMap = generateMap(nativePdb,nativeMtz)
        if not os.path.isfile(nativeMap):
            raise RuntimeError,"Cannot find nativeMap: {0}".format(nativeMap)
    else:
        shelxeExe=ample_util.find_exe('shelxe')
    
    for mrPdb in mrPdbs:
        if phenix:
            print "Searching for origin shift using: {0} {1}".format( nativeMap, mrPdb )
            origin =  phenixer.ccmtzOrigin( nativeMap, mrPdb )
            # offset.pdb is the mrPdb moved onto the new origin
            offsetPdb = "offset.pdb"
            print "Found origin: {0}\nOffset pdb is: {1}".format( origin, offsetPdb )
        else:
            originShift=shelxe.shelxe_origin(shelxeExe,nativePdb, nativeMtz, mrPdb)
            print "Found origin: {0}".format( originShift )
            offsetPdb=ample_util.filename_append(mrPdb, astr='offset', directory=os.getcwd())
            pdb_edit.translate(mrPdb, offsetPdb, originShift)
        
        # Run csymmatch to map offsetted pdb onto native
        csymmPdb = ample_util.filename_append( filename=mrPdb, astr="csymmatch", directory=outDir )
        print "Running csymmatch to wrap {0} onto native {1}".format( offsetPdb, nativePdb )
        csymmatch.Csymmatch().run( refPdb=nativePdb, inPdb=offsetPdb, outPdb=csymmPdb, originHand=False )
        
        print "Matched PDB is: {0}".format( csymmPdb )
    
    return

def generateMap(nativePdb,nativeMtz):
    
    # Get the labels from the MTZ file
    print "Parsing MTZ file {0} to determine column labels".format( nativeMtz )
    F,SIGF,FREE=mtz_util.get_labels(nativeMtz)
    # Generate map from 
    print "Generating map from: {0} {1}".format( nativeMtz,
                                                 nativePdb )
    return phenixer.generateMap(nativeMtz,
                                nativePdb,
                                FP=F,
                                SIGFP=SIGF,
                                FREE=FREE)
    
if __name__ == "__main__":
    assert len(sys.argv) >= 4,"Usage: {0} native.pdb native.mtz molecular_replacement.pdb[s]".format( sys.argv[0] )
    nativePdb = sys.argv[1]
    nativeMtz = sys.argv[2]
    mrPdbs    = sys.argv[3:]
    run(nativePdb=nativePdb, nativeMtz=nativeMtz, mrPdbs=mrPdbs)

