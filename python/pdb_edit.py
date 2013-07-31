'''
Useful manipulations on PDB files
'''

import copy
import re
import types
import unittest

three2one = {
    'ALA' : 'A',    
    'ARG' : 'R',    
    'ASN' : 'N',    
    'ASP' : 'D',    
    'CYS' : 'C',    
    'GLU' : 'E',    
    'GLN' : 'Q',    
    'GLY' : 'G',    
    'HIS' : 'H',    
    'ILE' : 'I',    
    'LEU' : 'L',    
    'LYS' : 'K',    
    'MET' : 'M',    
    'PHE' : 'F',    
    'PRO' : 'P',    
    'SER' : 'S',    
    'THR' : 'T',    
    'TRP' : 'W',    
    'TYR' : 'Y',   
    'VAL' : 'V,'
}

# http://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
#aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )
one2three =  dict((v, k) for (k, v) in three2one.items()) 

class PDBEdit(object):
    """Class for editing PDBs
    
    """
    
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
        pdb_in.close()
        
        return

    def keep_matching( self, refpdb, targetpdb, outpdb ):
        """Create a new pdb file that only contains that atoms in targetpdb that are
        also in refpdb
        
        NB: Assumes that both pdb files only contain one model and one chain
        
        Args:
        refpdb: path to pdb that contains the minimal set of atoms we want to keep
        targetpdb: path to the pdb that will be stripped of non-matching atoms
        outpdb: output path for the stripped pdb
        """
    
        
        def _write_matching_residues( ref_residues, target_residues, outfh ):
            
            # Loop over each residue in turn
            for idx, atoms_and_lines  in sorted( target_residues.items() ):
                
                # Get ordered list of the ref atom names for this residue
                rnames = [ x.name for x in ref_residues[ idx ] ]
                
                # Remove any not matching
                atoms = []
                atom_lines = []
                for i, a in enumerate( atoms_and_lines[0] ):
                    if a.name in rnames:
                        atoms.append( atoms_and_lines[0][i] )
                        atom_lines.append( atoms_and_lines[1][i] )
                
                
                # Now just have matching so output in the correct order
                for refname in rnames:
                    for i, atom in enumerate( atoms ):
                        if atom.name == refname:
                            # Found the matching atom so write out the corresponding line
                            outfh.write( atom_lines[i] )
                            # now delete both this atom and the line
                            atoms.pop(i)
                            atom_lines.pop(i)
                            # jump out of inner loop
                            break
            return
    
        # Go through refpdb and find which ref_residues are present
        f = open(refpdb, 'r')
        
        # map of resSeq to list of PdbAtom objects for the reference residues
        ref_residues = {}
        
        last = None
        for line in f:
            if line.startswith("MODEL"):
                raise RuntimeError, "Multi-model file!"
            
            if line.startswith("ATOM"):
                a = PdbAtom( line )
                if a.resSeq != last:
                    if a.resSeq in ref_residues:
                        raise RuntimeError,"Multiple chains in pdb - found residue #: {0} again.".format(a.resSeq)
                    last = a.resSeq
                    ref_residues[ last ] = [ a ]
                else:
                    ref_residues[ last ].append( a )
                    
        f.close()
        
        # Now read in target pdb and output everything bar the atoms in this file that
        # don't match those in the refpdb
        t = open(targetpdb,'r')
        out = open(outpdb,'w')
        
        reading=-1 # The residue we are reading - set to -1 when we are not reading
        
        target_residues = {} # dict mapping residue index to a a tuple of (atoms, lines), where atoms is a list of the atom
        # objects and lines is a list of the lines used to create the atom objects
        
        for line in t:
            
            if line.startswith("MODEL"):
                raise RuntimeError, "Multi-model file!"

            if line.startswith("ANISOU"):
                raise RuntimeError, "I cannot cope with ANISOU! {0}".format(line)
            
            # Stop at TER
            if line.startswith("TER"):
                # Write out the matching residues
                _write_matching_residues( ref_residues, target_residues, out )
                out.write("TER\n\n")
                break
            
            if line.startswith("ATOM"):
                
                atom = PdbAtom( line )
                
                # We copy resSeq to make sure we don't use a reference for our index
                resSeq = copy.copy( atom.resSeq )
                
                # Skip any ref_residues that don't match
                if resSeq in ref_residues:
                
                    # If this is the first one add the empty tuple and reset reading
                    if reading != resSeq:
                        # each tuple is a list of atom objects and lines
                        target_residues[ resSeq ] = ( [], [] )
                        reading = resSeq
                        
                    target_residues[ resSeq ][0].append( atom )
                    target_residues[ resSeq ][1].append( line )
                    
                # we don't write out any atom lines as they are either not matching or 
                # we write out matching at the end
                continue
            
            #Endif line.startswith("ATOM")
            
            # Output everything else
            out.write(line)
            
        t.close()
        out.close()
        
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
        
        return
    
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
    
    def strip_hetatm( self, inpath, outpath):
        """Remove all hetatoms from pdbfile"""
        o = open( outpath, 'w' )
        
        hremoved=-1
        for i, line in enumerate( open(inpath) ):
            
            # Remove EOL
            line = line.rstrip( "\n" )
            
            # Remove any HETATOM lines and following ANISOU lines
            if line.startswith("HETATM"):
                hremoved = i
                continue
            
            if line.startswith("ANISOU") and i == hremoved+1:
                continue
            
            o.write( line + "\n" )
            
        o.close()
        
        return
    
    def to__std_chain( self, pdbin, pdbout ):
        """Take the first chain of the first model and switch any non-standard AA's to their standard names.
        We also remove any ANISOU lines.
        """
    
        modres = [] # List of modres objects
        modres_names = []
        gotModel=False # to make sure we only take the first model
        chain=None # make sure we only get one chain
        reading=False # If reading structure
        
        
        pdbinf = open(pdbin,'r')
        pdboutf = open(pdbout,'w')
        line = True # Just for the first line
        while line:
    
            # Read in the line
            line = pdbinf.readline()
                    
            # Skip any ANISOU lines
            if line.startswith("ANISOU"):
                continue
            
            # Extract all MODRES DATA
            if line.startswith("MODRES"):
                modres.append( PdbModres( line ) )
                
            # Only extract the first model
            if line.startswith("MODEL"):
                if gotModel:
                    raise RuntimeError,"Found additional model! {0}".format( line )
                else:
                    gotModel=True
            
            # End at the end of the first chain
            if line.startswith("TER"):
                pdboutf.write( line )
                pdboutf.write( "END\n" )
                break
            
            # First time we hit coordinates we set up our data structures
            if not reading and ( line.startswith("HETATM") or line.startswith("ATOM") ):
                # There is a clever way to do this with list comprehensions but this is not it...
                for m in modres:
                    if m.resName not in modres_names:
                        modres_names.append( m.resName )
                        
                # Now we're reading
                reading=True
                    
            # Switch any residue names
            if len( modres):
                if line.startswith("HETATM"):
                    
                    hetatm = PdbHetatm( line )
                    
                    # See if this HETATM is in the chain we are reading and one of the residues to change
                    if hetatm.resName in modres_names:
                        for m in modres:
                            if hetatm.resName == m.resName and hetatm.chainID == m.chainID:
                                # Change this HETATM to an ATOM
                                atom = PdbAtom().fromHetatm( hetatm )
                                # Switch residue name
                                atom.resName = m.stdRes
                                # Convert to a line
                                line = atom.toLine()+"\n"
                                break
            
            # Any HETATM have been dealt with so just process as usual
            if line.startswith("ATOM"):
                atom = PdbAtom( line )
                if not chain:
                    chain = atom.chainID
                    
                if chain and atom.chainID != chain:
                    raise RuntimeError, "Got a different chain! {0}".format(line)
                
                if atom.resName not in three2one:
                    raise RuntimeError, "Unrecognised residue! {0}".format(line)
                    
            # Output everything else
            pdboutf.write( line )
    
            # END reading loop
            
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
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
"""
    def __init__(self, line=None):
        """Set up attributes"""
        
        if line:
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
        self.name = line[12:16].strip()
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
        if len(line) >= 77 and line[76:78].strip():
            self.element = line[76:78].strip()
        if len(line) >= 80 and line[78:80].strip():
            self.charge = int(line[78:80])
    
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
        self.name = line[12:16].strip()
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
        if len(line) >= 77 and line[76:78].strip():
            self.element = line[76:78].strip()
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
        self.assertEqual(a.name,'NH1')
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
        self.assertEqual(a.name,'MG')
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
    
    def testWriteHetatm(self):
        """Round-trip an atom line"""
        
        line = "HETATM 8239   O1 SO4 A2001      11.191 -14.833 -15.531  1.00 50.12           O  "
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
        self.assertEqual(a.stdRes,'G')
        self.assertEqual(a.comment,"1N-METHYLGUANOSINE-5'-MONOPHOSPHATE")
    
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
   
           
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
