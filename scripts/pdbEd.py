#!/usr/bin/env python
'''
Standalone scripts to manipulate PDB files
'''

# python imports
import argparse
import os
import sys
import tempfile


# Get the path to the ample python directory and add it to our path
apdir = os.path.abspath( os.path.join( os.path.dirname(__file__),"../python" ) )
sys.path.append( apdir )

# our imports
import ample_util
import pdb_edit


#import logging
#logging.basicConfig()
#logging.getLogger().setLevel(logging.DEBUG)

def extract_chain( inpdb, outpdb, chainID=None, newChainID=None ):
    """Extract chainID from inpdb and renumner"""
    
    if not newChainID:
        newChainID = chainID
    
    logfile = outpdb+".log"
    cmd="pdbcur xyzin {0} xyzout {1}".format( inpdb, outpdb ).split()
    stdin="""lvchain {0}
renchain {0} {1}
sernum
""".format( chainID, newChainID )
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
    
    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
        
    return

def get_info( inpdb ):
    PE = pdb_edit.PDBEdit()
    info = PE.get_info( inpdb )
    return info

def keep_matching( refpdb, targetpdb, outpdb ):
    """Only keep those atoms in targetpdb that are in refpdb and write the result to outpdb.
    We also take care of renaming any chains.
    """

    PE = pdb_edit.PDBEdit()

    # First get info on the two models
    refinfo = PE.get_info( refpdb )
    if len(refinfo.models) > 1:
        raise RuntimeError, "refpdb {0} has > 1 model!".format( refpdb )
    
    targetinfo = PE.get_info( targetpdb )
    if len(targetinfo.models) > 1:
        raise RuntimeError, "targetpdb {0} has > 1 model!".format( targetpdb )
    
    # If the chains have different names we need to rename the target to match the reference
    if refinfo.models[0].chains != targetinfo.models[0].chains:
        #print "keep_matching CHAINS ARE DIFFERENT BETWEEN MODELS: {0} : {1}".format(refinfo.models[0].chains, targetinfo.models[0].chains )
        
        if len(refinfo.models[0].chains) != len(targetinfo.models[0].chains):
            raise RuntimeError, "Different numbers of chains!"
        
        # We need to rename all the chains target to match those in refpdb using pdbcur
        targettmp = tmpFileName()+".pdb" # pdbcur insists names have a .pdb suffix
        
        stdint = ""
        for i, refchain in enumerate( refinfo.models[0].chains ):
            stdint += "renchain  /*/{0} '{1}'\n".format( targetinfo.models[0].chains[i], refchain )
 
        # now renumber with pdbcur
        logfile = targettmp+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( targetpdb, targettmp ).split()
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=True, stdin=stdint)
        
        if retcode == 0:
            # remove temporary files
            os.unlink(logfile)
        else:
            raise RuntimeError,"Error renaming chains!"   
        
    # Now we do our keep matching    
    tmp1 = tmpFileName()+".pdb" # pdbcur insists names have a .pdb suffix 
    
    PE = pdb_edit.PDBEdit()
    PE.keep_matching( refpdb, targettmp, tmp1 )
    
    # now renumber with pdbcur
    logfile = tmp1+".log"
    cmd="pdbcur xyzin {0} xyzout {1}".format( tmp1, outpdb ).split()
    stdint="""sernum
"""
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdint)
    
    if retcode == 0:
        # remove temporary files
        os.unlink(tmp1)
        os.unlink(targettmp)
        os.unlink(logfile)
    
    return retcode

def reforigin_rmsd( refpdb, inpdb, outpdb, DMAX=100 ):
    
    logfile = outpdb+".log"
    cmd="reforigin xyzin {0} xyzref {1} xyzout {2} DMAX {3}".format( refpdb, inpdb, outpdb, DMAX ).split()
    
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False)
    
    if retcode != 0:
        raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
    
    # Parse logfile to get RMSD
    rms = None
    for line in open( logfile, 'r' ):
        if line.startswith("RMS deviation:"):
            rms = float( line.split()[-1] )
            break
    
    if not rms:
        raise RuntimeError, "Error extracting RMS from logfile: {0}".format( logfile )
    
    return rms

def tmpFileName():
    """Return a filename for a temporary file"""

    # Get temporary filenames
    t = tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True)
    tmp1 = t.name
    t.close()
    return tmp1
    
def standardise( inpdb, outpdb, model=None):
    """Rename any non-standard AA, remove solvent and only keep most probably conformation.
    
    Args:
    model: extract this model from the file
    """
    

    tmp1 = tmpFileName()
    tmp1+=".pdb" # pdbcur insists names have a .pdb suffix
    
    # Strip to one chain and standardise AA names
    PE = pdb_edit.PDBEdit()
    PE.std_residues(inpdb, tmp1)
    
    # Now clean up with pdbcur
    logfile = tmp1+".log"
    cmd="pdbcur xyzin {0} xyzout {1}".format( tmp1, outpdb ).split()
    stdin="""delsolvent
mostprob
"""
    
    if model:
        assert type(model) == int
        stdin += "lvmodel /{0}\n".format( model )
        
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
    
    if retcode == 0:
        # remove temporary files
        os.unlink(tmp1)
        os.unlink(logfile)
    
    return retcode

#
# Command-line handling
#
parser = argparse.ArgumentParser(description='Manipulate PDB files', prefix_chars="-")

group = parser.add_mutually_exclusive_group()
group.add_argument('-one_std_chain', action='store_true',
                   help='Take pdb to one model/chain that contains only standard amino acids')

group.add_argument('-keep_matching', action='store_true',
                   help='keep matching atoms')

parser.add_argument('-ref_file', type=str,
                   help='The reference file')

parser.add_argument('input_file',
                   help='The input file - will not be altered')

parser.add_argument('output_file',
                   help='The output file - will be created')

#refpdb="/home/jmht/Documents/test/3PCV/test/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD_chain1.pdb"
#targetpdb="/home/jmht/Documents/test/3PCV/test/3PCV_clean.pdb"
#outpdb="/home/jmht/Documents/test/3PCV/test/matching1.pdb"

# refpdb="/home/jmht/Documents/test/3U2F/molrep/refine/refmac_molrep_loc0_ALL_poly_ala_trunc_0.21093_rad_2_UNMOD.pdb"
# targetpdb="/home/jmht/Documents/test/3U2F/test/3U2F_clean.pdb"
# outpdb="/home/jmht/Documents/test/3U2F/test/matching.pdb"
# keep_matching( refpdb, targetpdb, outpdb )

#to_1_std_chain("/home/jmht/Documents/test/3U2F/3U2F.pdb","/home/jmht/Documents/test/3U2F/test/3U2F_clean.pdb")


if "__name__" == "__main__":
    args = parser.parse_args()
    
    # Get full paths to all files
    args.input_file = os.path.abspath( args.input_file )
    if not os.path.isfile(args.input_file):
        raise RuntimeError, "Cannot find input file: {0}".format( args.input_file )
    args.output_file = os.path.abspath( args.output_file )
    if args.ref_file:
        args.ref_file = os.path.abspath( args.ref_file )
        if not os.path.isfile(args.ref_file):
            raise RuntimeError, "Cannot find ref file: {0}".format( args.ref_file )
    
    if args.one_std_chain:
        to_1_std_chain( args.input_file, args.output_file )
    elif args.keep_matching:
        keep_matching( args.ref_file, args.input_file, args.output_file )