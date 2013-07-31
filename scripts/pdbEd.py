#!/usr/bin/env python
'''
Standalone scripts to manipulate PDB files
'''

# python imports
import os
import tempfile

# our imports
import ample_util
import pdb_edit


#import logging
#logging.basicConfig()
#logging.getLogger().setLevel(logging.DEBUG)


def to_1_std_chain( inpdb, outpdb):
    """Clean down to 1 model/chain with standard amino acids"""
    
    
    # Get temporary filenames
    t = tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=False)
    tmp1 = t.name
    t.close()
    os.unlink(tmp1)
    tmp1+=".pdb" # pdbcur insists names have a .pdb suffix
    
    # Strip to one chain and standardise AA names
    PE = pdb_edit.PDBEdit()
    PE.to_1_std_chain(inpdb, tmp1)
    
    # Now clean up with pdbcur
    cmd="/opt/ccp4-6.3.0/bin/pdbcur xyzin {0} xyzout {1}".format( tmp1, outpdb ).split()
    stdin="""delsolvent
mostprob
"""
    retcode = ample_util.run_command(cmd=cmd, logfile=None, directory=os.getcwd(), dolog=False, stdin=stdin)
    
    # remove temporary file
    os.unlink(tmp1)
    
    return retcode

# root="/home/jmht/Documents/test/3PCV"
# os.chdir(root)
# to1stdChain(root+"/3PCV.pdb",root+"/test/3PCV_clean.pdb")

def keep_matching( refpdb, targetpdb, outpdb ):
    PE = pdb_edit.PDBEdit()
    PE.keep_matching( refpdb, targetpdb, outpdb )
    
    # now renumber

refpdb="/home/jmht/Documents/test/3PCV/test/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD_chain2.pdb"
targetpdb="/home/jmht/Documents/test/3PCV/test/3PCV_clean.pdb"
outpdb="/home/jmht/Documents/test/3PCV/test/matching2.pdb"
keep_matching( refpdb, targetpdb, outpdb )
