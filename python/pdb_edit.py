'''
Created on 15 Apr 2013

@author: jmht
'''

import re

class PDBEdit(object):
    """Class for editing PDBs
    
    """
    
    def select_residues(self, infile=None, outfile=None, residues=None ):
        """Create a new pdb by selecting only the numbered residues from the list.
        
        Args:
        infile: path to input pdb
        outfile: path to output pdb
        residues: list of integers of the residues to keep
        
        Return:
        path to new pdb or None
        """
    
        assert infile, outfile
        assert type(residues) == list
    
        pdb_in = open(infile, "r")
        pdb_out = open(outfile , "w")
        
        # Loop through PDB files and create new ones that only contain the residues specified in the list
        for pdbline in pdb_in:
            pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
            pdb_result = pdb_pattern.match(pdbline)
            if pdb_result:
                pdb_result2 = re.split(pdb_pattern, pdbline )
                for i in residues : #convert to ints to comparex
        
                    if int(pdb_result2[5]) == int(i):
        
                        pdb_out.write(pdbline)
        
        pdb_out.close()
        
        return outfile