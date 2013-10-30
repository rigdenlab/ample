'''
Created on 7 Aug 2013

@author: jmht

Classes for holding data from PDB files
'''

import copy
import os
import types
import unittest

# Non-redundant origins from:
# http://www.ccp4.ac.uk/dist/html/alternate_origins.html
_origins = {  
                   
                   # TRICLINIC
                   '1aP' : [ 
                             [ 0.0, 0.0, 0.0 ],
                             ],
                   
                   # MONOCLINIC
                   '2mP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           [ 0.5, 0.0, 0.0 ],
                           [ 0.5, 0.0, 0.5 ],
                           ] ,
  
                   '2mC' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           ] ,      
                    
                   '2mA' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.5, 0.0, 0.0 ],
                           ] ,
                    
                   '2mI' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           ] ,
                   
                   # ORTHORHOMBIC
                    '222oP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           [ 0.0, 0.5, 0.0 ],
                           [ 0.0, 0.5, 0.5 ],
                           [ 0.5, 0.0, 0.0 ],
                           [ 0.5, 0.0, 0.5 ],
                           [ 0.5, 0.5, 0.0 ],
                           [ 0.5, 0.5, 0.5 ],
                           ] ,
                                  
                    '222oC' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           [ 0.5, 0.0, 0.0 ],
                           [ 0.5, 0.0, 0.5 ],
                           ] ,     
                    
                    '222oF' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.25, 0.25, 0.25 ],
                           [ 0.5, 0.5, 0.5 ],
                           [ 0.75, 0.75, 0.75 ],
                           ] ,  
                    
                    '222oI' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           [ 0.0, 0.5, 0.0 ],
                           [ 0.5, 0.0, 0.0 ],
                           ] ,
                   
                   # TETRAGONAL
                    '4tP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.5, 0.5, 0.0 ],
                           ] ,           
                    
                    '4tI' : [
                           [ 0.0, 0.0, 0.0 ],
                           ] ,             

                    '422tP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           [ 0.5, 0.5, 0.0 ],
                           [ 0.5, 0.5, 0.5 ],
                           ] ,
                   
                    '422tI' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           ] ,
                   
                   # TETRAGONAL
                   '4tP' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.5, 0.5, 0.0 ],
                          ] ,
                   
                   '4tI' : [
                          [ 0.0, 0.0, 0.0 ],
                          ] ,
                                      
                   '422tP' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.0, 0.0, 0.5 ],
                          [ 0.5, 0.5, 0.0 ],
                          [ 0.5, 0.5, 0.5 ],
                          ] ,
                   
                   '422tI' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.0, 0.0, 0.5 ],
                          ] ,
                   
                   # TRIGONAL
                   '3hP' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ float(1/3), float(2/3), 0.0 ],
                          [ float(2/3), float(1/3), 0.0 ],
                          ] ,
                   
                   '3hR_1' : [
                          [ 0.0, 0.0, 0.0 ],
                          ] ,
                   
                   '3hR_2' : [
                          [ 0.0, 0.0, 0.0 ],
                          ] ,
                   
                   '312hP' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.0, 0.0, 0.5 ],
                          [ float(1/3), float(2/3), 0.0 ],
                          [ float(1/3), float(2/3), 0.5 ],
                          [ float(2/3), float(1/3), 0.0 ],
                          [ float(2/3), float(1/3), 0.5 ],
                          ] ,
                   
                   '321hP' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.0, 0.0, 0.5 ],
                          ] ,
                   
                   '32hR_1' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.0, 0.0, 0.5 ],
                          ] ,
                   
                   '32hR_2' : [
                          [ 0.0, 0.0, 0.0 ],
                          [ 0.5, 0.5, 0.5 ],
                          ] ,
                                      
                    # HEXAGONAL
                    '6hP' : [
                           [ 0.0, 0.0, 0.0 ],
                           ] ,
                          
                    '622hP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.0, 0.0, 0.5 ],
                           ] ,
                          
                    # CUBIC
                    '23cP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.5, 0.5, 0.5 ],
                           ] ,
                                          
                    '23cF' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.25, 0.25, 0.25 ],
                           [ 0.5, 0.5, 0.5 ],
                           [ 0.75, 0.75, 0.75 ],
                           ] ,

                    '23cI' : [
                           [ 0.0, 0.0, 0.0 ],
                           ] ,                                  

                    '432cP' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.5, 0.5, 0.5 ],
                           ] ,
                                              
                    '432cF' : [
                           [ 0.0, 0.0, 0.0 ],
                           [ 0.5, 0.5, 0.5 ],
                           ] ,                                  
                                              
                    '432cI' : [
                           [ 0.0, 0.0, 0.0 ],
                           ] ,                                  
                    }

_spacegroup2origin = {
              # Primitive
              'P1'          : _origins[ '1aP' ],
              # MONOCLINIC
              'P2'          : _origins[ '2mP' ],
              
              'P21'         : _origins[ '2mP' ],
              
              'C2'          : _origins[ '2mC' ],
              
              'A2'          : _origins[ '2mA' ],
              
              'I2'          : _origins[ '2mI' ],
              
              # ORTHORHOMBIC
              'P 2 2 2'      : _origins[ '222oP' ],
              'P 21 2 2'     : _origins[ '222oP' ],
              'P 2 21 2'     : _origins[ '222oP' ],
              'P 2 2 21'     : _origins[ '222oP' ],
              'P 2 21 21'    : _origins[ '222oP' ],
              'P 21 2 21'    : _origins[ '222oP' ],
              'P 21 21 2'    : _origins[ '222oP' ],
              'P 21 21 21'   : _origins[ '222oP' ],
              
              'C 2 2 21'     : _origins[ '222oC' ],
              'C 2 2 2'      : _origins[ '222oC' ],
              
              'F 2 2 2'      : _origins[ '222oF' ],
              
              'I 2 2 2'      : _origins[ '222oI' ],
              'I 21 12 21'   : _origins[ '222oI' ],
              
              # TETRAGONAL
              'P 4'          : _origins[ '4tP' ],
              'P 4 1'        : _origins[ '4tP' ],
              'P 4 2'        : _origins[ '4tP' ],
              'P 4 3'        : _origins[ '4tP' ],
              
              'I 4'          : _origins[ '4tI' ],
              'I 4 1'        : _origins[ '4tI' ],
              
              'P 4 2 2'      : _origins[ '422tP' ],
              'P 4 21 2'     : _origins[ '422tP' ],
              'P 41 2 2'     : _origins[ '422tP' ],
              'P 41 21 2'    : _origins[ '422tP' ],
              'P 42 2 2'     : _origins[ '422tP' ],
              'P 42 21 2'    : _origins[ '422tP' ],
              'P 43 2 2'     : _origins[ '422tP' ],
              'P 43 21 2'    : _origins[ '422tP' ],
              
              'I 4 2 2'      : _origins[ '422tI' ],
              'I 41 2 2'     : _origins[ '422tI' ],
              
              # TRIGONAL
              'P 3'          : _origins[ '3hP' ],
              'P 31'         : _origins[ '3hP' ],
              'P 32'         : _origins[ '3hP' ],
              
              'H 3'          : _origins[ '3hR_1' ],
              'R 3'          : _origins[ '3hR_2' ],
              
              'P 3 1 2'      : _origins[ '312hP' ],
              'P 31 1 2'     : _origins[ '312hP' ],
              'P 32 1 2'     : _origins[ '312hP' ],
              
              'P 3 2 1'     : _origins[ '321hP' ],
              'P 31 2 1'    : _origins[ '321hP' ],
              'P 32 2 1'    : _origins[ '321hP' ],
              
              'H 3 2'       : _origins[ '32hR_1' ],
              'R 3 2'       : _origins[ '32hR_2' ],
              
              # HEXAGONAL
              'P 6'         : _origins[ '6hP' ],
              'P 6 1'       : _origins[ '6hP' ],
              'P 6 5'       : _origins[ '6hP' ],
              'P 6 2'       : _origins[ '6hP' ],
              'P 6 4'       : _origins[ '6hP' ],
              'P 6 3'       : _origins[ '6hP' ],
              
              'P 6 2 2'     : _origins[ '622hP' ],
              'P 61 2 2'    : _origins[ '622hP' ],
              'P 65 2 2'    : _origins[ '622hP' ],
              'P 62 2 2'    : _origins[ '622hP' ],
              'P 64 2 2'    : _origins[ '622hP' ],
              'P 63 2 2'    : _origins[ '622hP' ],
              
              # CUBIC
              'P 2 3'       : _origins[ '23cP' ],
              'P 21 3'      : _origins[ '23cP' ],
              
              'F 2 3'       : _origins[ '23cF' ],
              
              'I 2 3'       : _origins[ '23cI' ],
              'I 21 3'      : _origins[ '23cI' ],
              
              'P 4 3 2'     : _origins[ '432cP' ],
              'P 42 3 2'    : _origins[ '432cP' ],
              'P 43 3 2'    : _origins[ '432cP' ],
              'P 41 3 2'    : _origins[ '432cP' ],
              
              'F 4 3 2'     : _origins[ '432cF' ],
              'F 41 3 2'    : _origins[ '432cF' ],
              
              'I 4 3 2'     : _origins[ '432cI' ],
              'I 41 3 2'    : _origins[ '432cI' ],
              
              }


#symoplib = "/Applications/ccp4-6.4.0/lib/data/symop.lib"
def _altlabel( spaceGroup, symoplib=None ):
    
    if not symoplib:
        symoplib = os.path.join( os.environ['CCP4'], "lib/data/symop.lib" )
        
    for line in open( symoplib, 'r' ):
        if "'" in line:
            # Assume first single-quote enclosed string is the one we want
            i = line.index( "'" )
            j = line.index("'", i+1 )
            sg = line[ i+1:j ]
            if spaceGroup == sg:
                return line.split()[ 3 ]

    raise KeyError, spaceGroup
    return 


def alternateOrigins( spaceGroupLabel ):
    """Given a space group label, return a list of (non-redundant) alternate
    origins as a list of float triples"""
    
    try:
        origins = _spacegroup2origin[ spaceGroupLabel ]
    except KeyError:
        label = _altlabel( spaceGroupLabel )
        origins = _spacegroup2origin[ label ]
    
    # Need to return a copy or if we manipulate the origins, we manipulate the copy
    # in the module
    return copy.copy( origins )


class CrystalInfo(object):
    def __init__(self, line=None):
        """foo"""
        
        self._reset()
        
        if line:
            self.fromLine( line )
        
        return
    
    def _reset( self ):
        
        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.spaceGroup = None
        self.z = None
        
        return
    
    def fromLine(self, line ):
        
        self.a = float( line[6:15] )
        self.b = float( line[15:24] )
        self.c = float( line[24:33] )
        self.alpha = float( line[33:40] )
        self.beta = float( line[40:47] )
        self.gamma = float( line[47:54] )
        self.spaceGroup = line[55:66].strip()
        self.z = int( line[66:70] )
        
        return

class PdbInfo(object):
    """A class to hold information extracted from a PDB file"""
    
    def __init__(self ):
        
        self.models = [] # List of PdbModel objects
        
        self.title = None # First line of the title
        self.resolution = None
        
        # http://www.wwpdb.org/documentation/format33/remarks1.html#REMARK%20280
        self.solventContent = None
        self.matthewsCoefficient = None
        
        self.crystalInfo = None
        
        return
    
class PdbModel(object):
    """A class to hold information on a single model in a PDB file"""
    
    def __init__(self ):
        
        self.serial = None
        self.chains = [] # Ordered list of chain IDs
        
        self.resSeqs = [] # Ordered list of list of resSeqs for each chain - matches order in self.chains
        self.sequences = [] # Ordered list of list of sequences for each chain - matches order in self.chains
        self.caMask = [] # Ordered list of list of booleans of residues with no CA atoms - matches order in self.chains
        self.bbMask = [] # Ordered list of list of boleans of residues with no backbone atoms - matches order in self.chains
        
        return

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
73 - 76        LString(4)    segID        Segment identifier, left-justified.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
"""
    def __init__(self, line=None):
        """Set up attributes"""
        
        self._setAtomType()
        
        if line:
            self.fromLine( line )
        
        return
    
    def _setAtomType(self):
        """This gets overridden in HETATM - otherwise everything the same"""
        self._atomType = "ATOM  "
        return
    
    def _reset(self):
        
        self.line = None # the line we were created from
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
        self.segID = None
        self.element = None
        self.charge = None
        
        return
    
    def _sanityCheck( self, line ):
        assert line[0:6] == self._atomType,"Line did not begin with an {0} record!: {1}".format( self._atomType, line )
        assert len(line) >= 54,"Line length was: {0}\n{1}".format(len(line),line)
        return
        
    def fromLine(self,line):
        """Initialise from the line from a PDB"""
        
        self._sanityCheck(line)
        self._reset()

        self.line = line
        
        self.serial = int(line[6:11])
        self.name = line[12:16]
        # Use for all so None means an empty field
        if line[16].strip():
            self.altLoc = line[16]
        self.resName = line[17:20].strip()
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
        if len(line) >= 76 and line[72:76].strip():
            self.segID = line[72:76].strip()
        if len(line) >= 77 and line[76:78].strip():
            self.element = line[76:78].strip()
        if len(line) >= 80 and line[78:80].strip():
            try:
                self.charge = int(line[78:80])
            except:
                raise RuntimeError, "Error getting charge ({0}) from line: {1}".format( line[78:80], line )
    
    def toLine(self):
        """Create a line suitable for printing to a PDB file"""
        
        s = self._atomType # 1-6
        s += "{0:5d}".format( self.serial ) # 7-11
        s += " " # 12 blank
        if len(self.name) != 4:
            raise RuntimeError,"Name must be 4 characters long!"
        s += "{0:4}".format( self.name ) # 13-16
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
        s += "{0:8.3F}".format( self.x ) #31-38
        s += "{0:8.3F}".format( self.y ) #39-46
        s += "{0:8.3F}".format( self.z ) #47-54
        if not self.occupancy: # 55-60
            s += "      "
        else:
            s += "{0:6.2F}".format( self.occupancy )
        if not self.tempFactor: # 61-66
            s += "      "
        else:
            s += "{0:6.2F}".format( self.tempFactor )
        s += "      " # 67-72 blank
        if not self.segID: # 73-76
            s += "    "
        else:
            s += "{0:>4}".format( self.segID )
        if not self.element: #77-78
            s += "  "
        else:
            s += "{0:>2}".format( self.element )
        if not self.charge: #79-80
            s += "  "
        else:
            s += "{0:2d}".format( self.charge )
            
        return s
    
    def fromHetatm( self, hetatm ):
        """Create Atom from Hetatm"""
        
        self.serial = hetatm.serial
        self.name = hetatm.name
        self.altLoc = hetatm.altLoc
        self.resName = hetatm.resName
        self.chainID = hetatm.chainID
        self.resSeq = hetatm.resSeq
        self.iCode = hetatm.iCode
        self.x = hetatm.x
        self.y = hetatm.y
        self.z = hetatm.z
        self.occupancy = hetatm.occupancy
        self.tempFactor = hetatm.tempFactor
        self.segID = hetatm.segID
        self.element = hetatm.element
        self.charge = hetatm.charge
        
        return self
        
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))

class PdbHetatm( PdbAtom ):
    """Identical to PdbAtom but just with a different _atomType"""
    
    def _setAtomType(self):
        self._atomType = "HETATM"
        return

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
        
        return
        
    def fromLine(self,line):
        """Initialise from the line from a PDB"""
        
        assert line[0:6] == "MODRES","Line did not begin with an MODRES record!: {0}".format(line)
        
        self._reset()
        
        self.idCode = line[7:11]
        self.resName = line[12:15].strip()
        # Use for all so None means an empty field
        if line[16].strip():
            self.chainID = line[16]
        self.seqNum = int(line[18:22])
        if line[22].strip():
            self.iCode = line[22]
        self.stdRes = line[24:27].strip()
        if line[29:70].strip():
            self.comment = line[29:70].strip()
            
        return
    
    def toLine(self):
        """Create a line suitable for printing to a PDB file"""
        
        s = "MODRES" # 1-6
        s += " " # 7 blank
        s += "{0:4}".format( self.idCode ) # 8-11
        s += " " # 12 blank
        s += "{0:>3}".format( self.resName ) # 13-15
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
        s += "{0:>3}".format( self.stdRes ) # 25-27
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
        self.assertEqual(a.element,'N')
        
        return
    
    def testReadAtom2(self):
        """Round-trip an atom line"""
        
        line = "ATOM     28  C   ALA A  12     -27.804  -2.987  10.849  1.00 11.75      AA-- C "
        a = PdbAtom( line )
        self.assertEqual(a.serial,28)
        self.assertEqual(a.name,' C  ')
        self.assertEqual(a.altLoc,None)
        self.assertEqual(a.resName,'ALA')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,12)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,-27.804)
        self.assertEqual(a.y,-2.987)
        self.assertEqual(a.z,10.849)
        self.assertEqual(a.occupancy,1.00)
        self.assertEqual(a.tempFactor,11.75)       
        self.assertEqual(a.segID,'AA--')
        self.assertEqual(a.element,'C')
        
        return
           
    def testWriteAtom1(self):
        """Round-trip an atom line"""
        
        line = "ATOM     41  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N  "
        a = PdbAtom( line )
        self.assertEqual( a.toLine(), line )
        
        return
           
    def testReadHetatm(self):
        """See if we can read a hetatom line"""

        line = "HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG  "
        a = PdbHetatm( line )
        self.assertEqual(a.serial,8237)
        self.assertEqual(a.name,'MG  ')
        self.assertEqual(a.altLoc,None)
        self.assertEqual(a.resName,'MG')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,1001)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,13.872)
        self.assertEqual(a.y,-2.555)
        self.assertEqual(a.z,-29.045)
        self.assertEqual(a.occupancy,1.00)
        self.assertEqual(a.tempFactor,27.36)
        self.assertEqual(a.element,'MG')
        
        return
    
    def testWriteHetatm(self):
        """Round-trip an atom line"""
        
        line = "HETATM 8239   O1 SO4 A2001      11.191 -14.833 -15.531  1.00 50.12           O  "
        a = PdbHetatm( line )
        self.assertEqual( a.toLine(), line )
        
        return
  
    def testWriteHetatm2(self):
        """Round-trip an atom line"""
        
        line = "HETATM    7  SD  FME A   1     -60.099  -1.874   3.446  1.00216.81           S  "
        a = PdbHetatm( line )
        self.assertEqual( a.toLine(), line )
        
        return
   
    def testReadModres(self):
        """See if we can read a modres line"""

        line = "MODRES 1IL2 1MG D 1937    G  1N-METHYLGUANOSINE-5'-MONOPHOSPHATE"
        a = PdbModres( line )
        self.assertEqual(a.idCode,"1IL2")
        self.assertEqual(a.resName,'1MG')
        self.assertEqual(a.chainID,'D')
        self.assertEqual(a.seqNum,1937)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.stdRes,'G')
        self.assertEqual(a.comment,"1N-METHYLGUANOSINE-5'-MONOPHOSPHATE")
        
        return
    
    def testWriteModres(self):
        """Round-trip a modres line"""
        
        line = "MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE"
        a = PdbModres( line )
        self.assertEqual(a.idCode,"2R0L")
        self.assertEqual(a.resName,'ASN')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.seqNum,74)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.stdRes,'ASN')
        self.assertEqual(a.comment,"GLYCOSYLATION SITE")
        self.assertEqual( a.toLine(), line )
        
        return
    
    def testReadCrystalInfo(self):
        """See if we can read a cryst1 line"""

        line = "CRYST1  117.000   15.000   39.000  90.00  90.00  90.00 P 21 21 21    8 "
        a = CrystalInfo( line )
        self.assertEqual(a.a,117.000)
        self.assertEqual(a.b,15.000)
        self.assertEqual(a.c,39.000)
        self.assertEqual(a.alpha,90)
        self.assertEqual(a.beta,90)
        self.assertEqual(a.gamma,90)
        self.assertEqual(a.spaceGroup,"P 21 21 21")
        self.assertEqual(a.z,8)
        
        return
           
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

