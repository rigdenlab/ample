'''
Created on 15 Apr 2013

@author: jmht
'''

import re
import types
import unittest

class PDBEdit(object):
    """Class for editing PDBs
    
    """
    
    def select_residues(self, inpath=None, outpath=None, residues=None ):
        """Create a new pdb by selecting only the numbered residues from the list.
        
        Args:
        infile: path to input pdb
        outfile: path to output pdb
        residues: list of integers of the residues to keep
        
        Return:
        path to new pdb or None
        """
    
        assert inpath, outpath
        assert type(residues) == list
    
        pdb_in = open(inpath, "r")
        pdb_out = open(outpath , "w")
        
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
        
        return
    
    def reliable_sidechains(self, inpath=None, outpath=None ):
        """Only output non-backbone atoms for residues in the res_names list.
        """
        
        # Remove sidechains that are in res_names where the atom name is not in atom_names
        res_names = [ 'MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
        atom_names = [ 'N', 'CA', 'C', 'O', 'CB' ]

        #   print 'Found ',each_file
        pdb_in = open( inpath, "r" )
        pdb_out = open( outpath, "w" )
        
        for pdbline in pdb_in:
            pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
            pdb_result = pdb_pattern.match(pdbline)
            
            # Check ATOM line and for residues in res_name, skip any that are not in atom names
            if pdb_result:
                pdb_result2 = re.split(pdb_pattern, pdbline)
                if pdb_result2[3] in res_names and not pdb_result2[2] in atom_names:
                    continue
            
            # Write out everything else
            pdb_out.write(pdbline)
        
        #End for
        pdb_out.close()
        pdb_in.close()
        
    def backbone(self, inpath=None, outpath=None ):
        """Only output backbone atoms.
        """        
        
        atom_names = [ 'N', 'CA', 'C', 'O', 'CB' ]

        #   print 'Found ',each_file
        pdb_in = open( inpath, "r" )
        pdb_out = open( outpath, "w" )    

        for pdbline in pdb_in:
            pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
            pdb_result = pdb_pattern.match(pdbline)
    
            if pdb_result:
                pdb_result2 = re.split(pdb_pattern, pdbline)
                if pdb_result2[3] != '':
                    if pdb_result2[2] not in atom_names:
                        continue
            
            # Write out everything else
            pdb_out.write(pdbline)
        
        #End for
        pdb_out.close()
        pdb_in.close



class PdbAtom(object):
    """
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
"""
    def __init__(self, line):
        """Set up attributes"""
        
        self.fromLine( line )
        
    
    def _reset(self):
        
        self.serial = None
        self.name = None
        self.altLoc = None
        self.resName = None
        self.chainID = None
        self.resSeq = None
        self.iCode = None
        self.x = None
        self.y = None
        self.z = None
        self.occupancy = None
        self.tempFactor = None
        self.element = None
        self.charge = None
        
    def fromLine(self,line):
        """Initialise from the line from a PDB"""
        
        
        assert line[0:6] == "ATOM  ","Line did not begin with an ATOM record!: {0}".format(line)
        assert len(line) >= 54,"Line length was: {0}\n{1}".format(len(line),line)
        
        self._reset()
        
        self.serial = int(line[6:11])
        self.name = line[12:16]
        # Use for all so None means an empty field
        if line[16].strip():
            self.altLoc = line[16]
        self.resName = line[17:20]
        if line[21].strip():
            self.chainID = line[21]
        if line[22:26].strip():
            self.resSeq = int(line[22:26])
        if line[26].strip():
            self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        if len(line) >= 60 and line[54:60].strip():
            self.occupancy = float(line[54:60])
        if len(line) >= 66 and line[60:66].strip():
            self.tempFactor = float(line[60:66])
        if len(line) >= 78 and line[76:78].strip():
            self.element = line[76:78]
        if len(line) >= 80 and line[78:80].strip():
            self.charge = line[78:80]
    
    def toLine(self):
        """Create a line suitable for printing to a PDB file"""
        
        s = "ATOM  " # 1-6
        s += "{0:5d}".format( self.serial ) # 7-11
        s += " " # 12 blank
        s += "{0:>4}".format( self.name ) # 13-16
        if not self.altLoc: #17
            s += " "
        else:
            s += "{0:1}".format( self.altLoc )
        s += "{0:3}".format( self.resName ) # 18-20
        s += " " # 21 blank
        if not self.chainID: #22
            s += " "
        else:
            s += "{0:1}".format( self.chainID )
        s += "{0:4}".format( self.resSeq ) #23-26
        if not self.iCode: #27
            s += " "
        else:
            s += "{0:1}".format( self.iCode )
        s += "   " # 28-30 blank
        s += "{0: 8.3F}".format( self.x ) #31-38
        s += "{0: 8.3F}".format( self.y ) #39-46
        s += "{0: 8.3F}".format( self.z ) #47-54
        if not self.occupancy: # 55-60
            s += "      "
        else:
            s += "{0: 6.2F}".format( self.occupancy )
        if not self.tempFactor: # 61-66
            s += "      "
        else:
            s += "{0: 6.2F}".format( self.tempFactor )
        s += "          " # 67-76 blank
        if not self.element: #77-78
            s += "  "
        else:
            s += "{0:>2}".format( self.element )
        if not self.charge: #79-80
            s += "  "
        else:
            s += "{0:2d}".format( self.charge )
            
        return s
        
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))



class PdbHetatm(object):
    """
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "HETATM"
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
"""
    def __init__(self, line):
        """Set up attributes"""
        
        self.fromLine( line )
        
    
    def _reset(self):
        
        self.serial = None
        self.name = None
        self.altLoc = None
        self.resName = None
        self.chainID = None
        self.resSeq = None
        self.iCode = None
        self.x = None
        self.y = None
        self.z = None
        self.occupancy = None
        self.tempFactor = None
        self.element = None
        self.charge = None
        
    def fromLine(self,line):
        """Initialise from the line from a PDB"""
        
        assert line[0:6] == "HETATM","Line did not begin with an HETATM record!: {0}".format(line)
        assert len(line) >= 54,"Line length was: {0}\n{1}".format(len(line),line)
        
        self._reset()
        
        self.serial = int(line[6:11])
        self.name = line[12:16]
        # Use for all so None means an empty field
        if line[16].strip():
            self.altLoc = line[16]
        self.resName = line[17:20]
        if line[21].strip():
            self.chainID = line[21]
        if line[22:26].strip():
            self.resSeq = int(line[22:26])
        if line[26].strip():
            self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        if len(line) >= 60 and line[54:60].strip():
            self.occupancy = float(line[54:60])
        if len(line) >= 66 and line[60:66].strip():
            self.tempFactor = float(line[60:66])
        if len(line) >= 78 and line[76:78].strip():
            self.element = line[76:78]
        if len(line) >= 80 and line[78:80].strip():
            self.charge = line[78:80]
    
    def toLine(self):
        """Create a line suitable for printing to a PDB file"""
        
        s = "HETATM" # 1-6
        s += "{0:5d}".format( self.serial ) # 7-11
        s += " " # 12 blank
        s += "{0:>4}".format( self.name ) # 13-16
        if not self.altLoc: #17
            s += " "
        else:
            s += "{0:1}".format( self.altLoc )
        s += "{0:3}".format( self.resName ) # 18-20
        s += " " # 21 blank
        if not self.chainID: #22
            s += " "
        else:
            s += "{0:1}".format( self.chainID )
        s += "{0:4}".format( self.resSeq ) #23-26
        if not self.iCode: #27
            s += " "
        else:
            s += "{0:1}".format( self.iCode )
        s += "   " # 28-30 blank
        s += "{0: 8.3F}".format( self.x ) #31-38
        s += "{0: 8.3F}".format( self.y ) #39-46
        s += "{0: 8.3F}".format( self.z ) #47-54
        if not self.occupancy: # 55-60
            s += "      "
        else:
            s += "{0: 6.2F}".format( self.occupancy )
        if not self.tempFactor: # 61-66
            s += "      "
        else:
            s += "{0: 6.2F}".format( self.tempFactor )
        s += "          " # 67-76 blank
        if not self.element: #77-78
            s += "  "
        else:
            s += "{0:>2}".format( self.element )
        if not self.charge: #79-80
            s += "  "
        else:
            s += "{0:2d}".format( self.charge )
            
        return s
        
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))


class PdbModres(object):
    """
COLUMNS        DATA TYPE     FIELD       DEFINITION
--------------------------------------------------------------------------------
 1 -  6        Record name   "MODRES"
 8 - 11        IDcode        idCode      ID code of this entry.
13 - 15        Residue name  resName     Residue name used in this entry.
17             Character     chainID     Chain identifier.
19 - 22        Integer       seqNum      Sequence number.
23             AChar         iCode       Insertion code.
25 - 27        Residue name  stdRes      Standard residue name.
30 - 70        String        comment     Description of the residue modification.
"""
    def __init__(self, line):
        """Set up attributes"""
        
        self.fromLine( line )
        
    
    def _reset(self):
        
        self.idCode = None
        self.resName = None
        self.chainID = None
        self.seqNum = None
        self.iCode = None
        self.stdRes = None
        self.comment = None
        
    def fromLine(self,line):
        """Initialise from the line from a PDB"""
        
        assert line[0:6] == "MODRES","Line did not begin with an MODRES record!: {0}".format(line)
        
        self._reset()
        
        self.idCode = line[7:11]
        self.resName = line[12:15]
        # Use for all so None means an empty field
        if line[16].strip():
            self.chainID = line[16]
        self.seqNum = int(line[18:22])
        if line[22].strip():
            self.iCode = line[22]
        self.stdRes = line[24:27]
        if line[29:70].strip():
            self.comment = line[29:70]
    
    def toLine(self):
        """Create a line suitable for printing to a PDB file"""
        
        s = "MODRES" # 1-6
        s += " " # 7 blank
        s += "{0:4}".format( self.idCode ) # 8-11
        s += " " # 12 blank
        s += "{0:3}".format( self.resName ) # 13-15
        s += " " # 16 blank
        if not self.chainID: #17
            s += " "
        else:
            s += "{0:1}".format( self.chainID )
        s += " " # 18 blank
        s += "{0:4d}".format( self.seqNum ) # 19-22
        if not self.iCode: #23
            s += " "
        else:
            s += "{0:1}".format( self.iCode )
        s += " " # 24 blank
        s += "{0:3}".format( self.stdRes ) # 25-27
        s += "  " # 28-29 blank
        if self.comment: # 30-70
            s += "{:<}".format( self.comment )
            
        return s
        
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))

class Test(unittest.TestCase):

    def testReadAtom(self):
        """See if we can read an atom line"""

        line = "ATOM     41  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N  "
        a = PdbAtom( line )
        self.assertEqual(a.serial,41)
        self.assertEqual(a.name,' NH1')
        self.assertEqual(a.altLoc,'A')
        self.assertEqual(a.resName,'ARG')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,-3)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,12.218)
        self.assertEqual(a.y,84.840)
        self.assertEqual(a.z,88.007)
        self.assertEqual(a.occupancy,0.5)
        self.assertEqual(a.tempFactor,40.76)
        self.assertEqual(a.element,' N')
    
    def testWriteAtom(self):
        """Round-trip an atom line"""
        
        line = "ATOM     41  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N  "
        a = PdbAtom( line )
        self.assertEqual( a.toLine(), line )
           
    def testReadHetatm(self):
        """See if we can read a hetatom line"""

        line = "HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG  "
        a = PdbHetatm( line )
        self.assertEqual(a.serial,8237)
        self.assertEqual(a.name,'MG  ')
        self.assertEqual(a.altLoc,None)
        self.assertEqual(a.resName,' MG')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,1001)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,13.872)
        self.assertEqual(a.y,-2.555)
        self.assertEqual(a.z,-29.045)
        self.assertEqual(a.occupancy,1.00)
        self.assertEqual(a.tempFactor,27.36)
        self.assertEqual(a.element,'MG')
    
    def testWriteHetatm(self):
        """Round-trip an atom line"""
        
        line = "HETATM 8239  O1  SO4 A2001      11.191 -14.833 -15.531  1.00 50.12           O  "
        a = PdbHetatm( line )
        self.assertEqual( a.toLine(), line )
   
    def testReadModres(self):
        """See if we can read a modres line"""

        line = "MODRES 1IL2 1MG D 1937    G  1N-METHYLGUANOSINE-5'-MONOPHOSPHATE"
        a = PdbModres( line )
        self.assertEqual(a.idCode,"1IL2")
        self.assertEqual(a.resName,'1MG')
        self.assertEqual(a.chainID,'D')
        self.assertEqual(a.seqNum,1937)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.stdRes,'  G')
        self.assertEqual(a.comment,"1N-METHYLGUANOSINE-5'-MONOPHOSPHATE")
    
    def testWriteModres(self):
        """Round-trip a modres line"""
        
        line = "MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE"
        a = PdbModres( line )
        self.assertEqual( a.toLine(), line )
   
           
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
