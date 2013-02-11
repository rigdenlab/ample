#!/usr/bin/env python
"""
edit the sidechains to make polyala, all and reliable
"""

import re
import sys
import shutil
import os

import unittest
import filecmp

########################################## ADD ALL

def edit_sidechains(each_file, outpath):
    """
    Prune down the ensembles based on their sidechain atoms
    All_atom: just copies the file
    SCWRL_Reliable_sidechains: for 'MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER'
      only output atom types: 'N', 'CA', 'C', 'O', 'CB'
    poly_ala: for all residues, just output the 'N', 'CA', 'C', 'O', 'CB'
    """

    if not os.path.exists(each_file):
        print 'Could not create ensemble for files (models too diverse): {}'.format(each_file)
        return
    
    #   print 'Found ',each_file
    my_infile = open (each_file)
    
    name = re.split('/', each_file)
    pdbname = str(name.pop()) 
    
    # For all atom just copy the file
    shutil.copy2( each_file,  outpath + 'All_atom_' + pdbname )
    
    # Remove sidechains that are in res_names where the atom name is not in atom_names
    res_names = [ 'MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
    atom_names = [ 'N', 'CA', 'C', 'O', 'CB' ]
    
    my_outfile = open (outpath + 'SCWRL_Reliable_sidechains_' + pdbname, "w")
    for pdbline in my_infile:
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)
        
        # Check ATOM line and for residues in res_name, skip any that are not in atom names
        if pdb_result:
            pdb_result2 = re.split(pdb_pattern, pdbline)
            if pdb_result2[3] in res_names and not pdb_result2[2] in atom_names:
                continue
        
        # Write out everything else
        my_outfile.write(pdbline)
    
    #End for
    my_outfile.close()

    # Rewind input file    
    my_infile.seek(0)
    my_outfile2 = open (outpath + 'poly_ala_' + pdbname, "w")
    # Only output mainchain atoms
    for pdbline in my_infile:
        #print pdbline
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)
        
        if pdb_result:
            pdb_result2 = re.split(pdb_pattern, pdbline)
            if pdb_result2[3] != '':
                if pdb_result2[2] not in atom_names:
                    continue
                
        my_outfile2.write(pdbline)
    
    my_outfile2.close()

class TestCell(unittest.TestCase):
    """
    Unit test
    """
    
    def testEditSidechains(self):
        """test we can prune sidechains"""
        
        root = "/opt/ample-dev1/"
        testdir = root + "tests/testfiles/"
        ensemble = testdir + "trunc_28.146439_rad_3.pdb"
        edit_sidechains(ensemble, testdir)
        
        allatom = "All_atom_trunc_28.146439_rad_3.pdb"
        reliable = "SCWRL_Reliable_sidechains_trunc_28.146439_rad_3.pdb"
        poly = "poly_ala_trunc_28.146439_rad_3.pdb"
        self.assertTrue( filecmp.cmp(testdir + allatom, 
                                     testdir + "trunc_28.146439_rad_3.pdb"),
                         "Error with All_atom" )
        self.assertTrue( filecmp.cmp(testdir + "orig." + reliable,
                                     testdir + reliable),
                         "Error with Reliable sidechains" )
        
        self.assertTrue( filecmp.cmp(testdir + "orig."+ poly,
                                     testdir + poly ),
                         "Error with polyalanine" )
        
        # Clean up
        os.unlink( testdir + allatom )
        os.unlink( testdir + reliable )
        os.unlink( testdir + poly )
        
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
