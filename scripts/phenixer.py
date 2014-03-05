'''
Created on 5 Mar 2014

@author: jmht
'''

import os
import sys
sys.path.append("/opt/ample-dev1/python")

import ample_util

def generateMap( mtz, pdb ):
    """Generate a map from an mtz file and a pdb using reforigin"""
    
    assert os.path.isfile( mtz ) and os.path.isfile( pdb ), "Cannot find files: {0} {1}".format( mtz, pdb )
    
    mapFile = ample_util.filename_append( filename=mtz, astr="map" )
    mapFile = os.path.abspath(mapFile)
    mapPdb = ample_util.filename_append( filename=pdb, astr="map" )

    cmd    = [ "refmac5", "HKLIN", mtz, "HKLOUT", mapFile, "XYZIN", pdb, "XYZOUT", mapPdb ]
    # FIX FOR DIFFERENT FP etc.     
    stdin ="""RIDG DIST SIGM 0.02
LABIN FP=FP SIGFP=SIGFP FREE=FREE
MAKE HYDR N
WEIGHT MATRIX 0.01
NCYC 0
END
"""
    ret = ample_util.run_command(cmd=cmd, logfile="generateMap.log", dolog=False, stdin=stdin)
    
    assert ret == 0, "generateMap refmac failed!"

    return mapFile

def ccmtzOrigin( nativeMap, mrPdb  ):
    """Use the phenix get_cc_mtz_pdb script to determine the origin of a MR pdb using the supplied map"""
    
    cmd = [ "phenix.get_cc_mtz_pdb", nativeMap, mrPdb ]
    ret = ample_util.run_command(cmd=cmd, logfile="get_cc_mtz_pdb.log", dolog=False )
    assert ret == 0, "phenix.get_cc_mtz_pdb refmac failed!"
    
    ofile = "temp_dir/resolve.offset"
    with open( ofile ) as o:
        line = o.readline().strip()
    
    t = line.split()
    assert t[0] == "OFFSET"
    origin = [ float( t[1] ) * -1, float( t[2] ) * -1, float( t[3] ) * -1 ]
    
    return origin


if __name__ == '__main__':
    
    os.chdir("/home/jmht/Documents/work/CC/phenix_get_cc_mtz_pdb/foo")
    mtz = "3T97-cad.mtz"
    native = "3T97_std.pdb"
    mrPdb = "phaser_loc0_ALL_poly_ala_trunc_2.733913_rad_1_UNMOD.1.pdb"
    mtzMap = generateMap( mtz, native )
    
    print ccmtzOrigin( mtzMap, mrPdb )