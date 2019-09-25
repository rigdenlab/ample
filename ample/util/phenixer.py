'''
Created on 5 Mar 2014

@author: jmht
'''

import os
import shutil
from ample.util import ample_util


def generateMap(mtz, pdb, FP='FP', SIGFP='SIGFP', FREE='FREE', directory=None):
    """Generate a map from an mtz file and a pdb using reforigin"""

    assert os.path.isfile(mtz) and os.path.isfile(pdb), "Cannot find files: {0} {1}".format(mtz, pdb)

    if not directory:
        directory = os.getcwd()

    mapFile = ample_util.filename_append(filename=mtz, astr="map", directory=directory)
    mapFile = os.path.abspath(mapFile)
    mapPdb = ample_util.filename_append(filename=pdb, astr="map", directory=directory)

    cmd = ["refmac5", "HKLIN", mtz, "HKLOUT", mapFile, "XYZIN", pdb, "XYZOUT", mapPdb]
    # FIX FOR DIFFERENT FP etc.
    stdin = """RIDG DIST SIGM 0.02
LABIN FP={0} SIGFP={1} FREE={2}
MAKE HYDR N
WEIGHT MATRIX 0.01
NCYC 0
END
""".format(
        FP, SIGFP, FREE
    )
    logfile = os.path.join(directory, "generateMap.log")
    ret = ample_util.run_command(cmd=cmd, logfile=logfile, dolog=True, stdin=stdin)

    assert ret == 0, "generateMap refmac failed-check log: {0}".format(logfile)

    return mapFile


def ccmtzOrigin(nativeMap, mrPdb):
    """Use the phenix get_cc_mtz_pdb script to determine the origin of a MR pdb using the supplied map"""

    # resolve can only handle file names < 75 characters so we need to truncate
    # We copy the file rather than symlink so that this works on windows and then delete afterwards
    tempnam = None
    if len(os.path.basename(mrPdb)) >= 75:
        tempnam = os.tempnam()
        # Need to add .pdb extension or it doesn't work
        tempnam += ".pdb"
        assert len(tempnam) < 75
        shutil.copy(mrPdb, tempnam)
        mrPdb = tempnam

    # make sure we can find the program
    get_cc_mtz_pdb = ample_util.find_exe('phenix.get_cc_mtz_pdb')
    cmd = [get_cc_mtz_pdb, nativeMap, mrPdb]
    ret = ample_util.run_command(cmd=cmd, logfile="get_cc_mtz_pdb.log", dolog=False)
    assert ret == 0, "phenix.get_cc_mtz_pdb refmac failed!"

    ofile = "temp_dir/resolve.offset"
    with open(ofile) as o:
        line = o.readline().strip()

    t = line.split()
    assert t[0] == "OFFSET"
    origin = [float(t[1]) * -1, float(t[2]) * -1, float(t[3]) * -1]

    # remove temp file if we created it
    if tempnam:
        os.unlink(tempnam)

    return origin


if __name__ == '__main__':

    os.chdir("/home/jmht/Documents/work/CC/phenix_get_cc_mtz_pdb/foo")
    mtz = "3T97-cad.mtz"
    native = "3T97_std.pdb"
    mrPdb = "phaser_loc0_ALL_poly_ala_trunc_2.733913_rad_1_UNMOD.1.pdb"
    mtzMap = generateMap(mtz, native)

    print (ccmtzOrigin(mtzMap, mrPdb))
