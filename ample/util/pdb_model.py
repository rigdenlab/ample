"""
Created on 7 Aug 2013

@author: jmht

Classes for holding data from PDB files
"""

import copy
import os
import types


class OriginInfo(object):
    def __init__(self, spaceGroupLabel=None):

        # These are reset on each call
        self._spaceGroup = None
        self._redundantSet = None
        self._nonRedundantSet = None
        self._floating = False

        self._setData()

        if spaceGroupLabel:
            self._getAlternateOrigins(spaceGroupLabel)

        return

    def _setData(self):
        # Non-redundant origins from:
        # http://www.ccp4.ac.uk/dist/html/alternate_origins.html
        # Organised in tuples, with True for the second item if the origin is one of the non-redundant set
        self._origins = {
            # TRICLINIC
            '1aP': [(['x', 'y', 'z'], True)],
            # MONOCLINIC
            '2mP': [([0.0, 'y', 0.0], True), ([0.0, 'y', 0.5], True), ([0.5, 'y', 0.0], True), ([0.5, 'y', 0.5], True)],
            '2mC': [
                ([0.0, 'y', 0.0], True),
                ([0.0, 'y', 0.5], True),
                ([0.5, 'y', 0.0], False),
                ([0.5, 'y', 0.5], False),
            ],
            '2mA': [
                ([0.0, 'y', 0.0], True),
                ([0.0, 'y', 0.5], False),
                ([0.5, 'y', 0.0], True),
                ([0.5, 'y', 0.5], False),
            ],
            '2mI': [
                ([0.0, 'y', 0.0], True),
                ([0.0, 'y', 0.5], True),
                ([0.5, 'y', 0.0], False),
                ([0.5, 'y', 0.5], False),
            ],
            # ORTHORHOMBIC
            '222oP': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([0.0, 0.5, 0.0], True),
                ([0.0, 0.5, 0.5], True),
                ([0.5, 0.0, 0.0], True),
                ([0.5, 0.0, 0.5], True),
                ([0.5, 0.5, 0.0], True),
                ([0.5, 0.5, 0.5], True),
            ],
            '222oC': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.5, 0.0], False),
                ([0.0, 0.5, 0.5], False),
                ([0.0, 0.0, 0.5], True),
                ([0.5, 0.0, 0.0], True),
                ([0.5, 0.0, 0.5], True),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], False),
            ],
            '222oF': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], False),
                ([0.0, 0.5, 0.0], False),
                ([0.0, 0.5, 0.5], False),
                ([0.25, 0.25, 0.25], True),
                ([0.25, 0.25, 0.75], False),
                ([0.25, 0.75, 0.25], False),
                ([0.25, 0.75, 0.75], False),
                ([0.5, 0.0, 0.0], False),
                ([0.5, 0.0, 0.5], False),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], True),
                ([0.75, 0.25, 0.25], False),
                ([0.75, 0.25, 0.75], False),
                ([0.75, 0.75, 0.25], False),
                ([0.75, 0.75, 0.75], True),
            ],
            '222oI': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([0.0, 0.5, 0.0], True),
                ([0.0, 0.5, 0.5], False),
                ([0.5, 0.0, 0.0], True),
                ([0.5, 0.0, 0.5], False),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], False),
            ],
            # TETRAGONAL
            '4tP': [([0.0, 0.0, 'z'], True), ([0.5, 0.5, 'z'], True)],
            '4tI': [([0.0, 0.0, 'z'], True), ([0.5, 0.5, 'z'], False)],
            '422tP': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([0.5, 0.5, 0.0], True),
                ([0.5, 0.5, 0.5], True),
            ],
            '422tI': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], False),
            ],
            # TRIGONAL
            '3hP': [
                ([0.0, 0.0, 'z'], True),
                ([float(1 / 3), float(2 / 3), 'z'], True),
                ([float(2 / 3), float(1 / 3), 'z'], True),
            ],
            '3hR_1': [
                ([0.0, 0.0, 'z'], True),
                ([float(1 / 3), float(2 / 3), 'z'], False),
                ([float(2 / 3), float(1 / 3), 'z'], False),
            ],
            '3hR_2': [(['x', 'x', 'x'], True)],
            '312hP': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([float(1 / 3), float(2 / 3), 0.0], True),
                ([float(1 / 3), float(2 / 3), 0.5], True),
                ([float(2 / 3), float(1 / 3), 0.0], True),
                ([float(2 / 3), float(1 / 3), 0.5], True),
            ],
            '321hP': [([0.0, 0.0, 0.0], True), ([0.0, 0.0, 0.5], True)],
            '32hR_1': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], True),
                ([float(1 / 3), float(2 / 3), float(1 / 6)], False),
                ([float(1 / 3), float(2 / 3), float(2 / 3)], False),
                ([float(2 / 3), float(1 / 3), float(1 / 3)], False),
                ([float(2 / 3), float(1 / 3), float(5 / 6)], False),
            ],
            '32hR_2': [([0.0, 0.0, 0.0], True), ([0.5, 0.5, 0.5], True)],
            # HEXAGONAL
            '6hP': [([0.0, 0.0, 'z'], True)],
            '622hP': [([0.0, 0.0, 0.0], True), ([0.0, 0.0, 0.5], True)],
            # CUBIC
            '23cP': [([0.0, 0.0, 0.0], True), ([0.5, 0.5, 0.5], True)],
            '23cF': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], False),
                ([0.0, 0.5, 0.0], False),
                ([0.0, 0.5, 0.5], False),
                ([0.25, 0.25, 0.25], True),
                ([0.25, 0.25, 0.75], False),
                ([0.25, 0.75, 0.25], False),
                ([0.25, 0.75, 0.75], False),
                ([0.5, 0.0, 0.0], False),
                ([0.5, 0.0, 0.5], False),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], True),
                ([0.75, 0.25, 0.25], False),
                ([0.75, 0.25, 0.75], False),
                ([0.75, 0.75, 0.25], False),
                ([0.75, 0.75, 0.75], True),
            ],
            '23cI': [([0.0, 0.0, 0.0], True), ([0.5, 0.5, 0.5], False)],
            '432cP': [([0.0, 0.0, 0.0], True), ([0.5, 0.5, 0.5], True)],
            '432cF': [
                ([0.0, 0.0, 0.0], True),
                ([0.0, 0.0, 0.5], False),
                ([0.0, 0.5, 0.0], False),
                ([0.0, 0.5, 0.5], False),
                ([0.5, 0.0, 0.0], False),
                ([0.5, 0.0, 0.5], False),
                ([0.5, 0.5, 0.0], False),
                ([0.5, 0.5, 0.5], True),
            ],
            '432cI': [([0.0, 0.0, 0.0], True), ([0.5, 0.5, 0.5], False)],
        }

        self._spacegroup2origin = {
            # Primitive
            'P1': self._origins['1aP'],
            # MONOCLINIC
            'P2': self._origins['2mP'],
            'P21': self._origins['2mP'],
            'C2': self._origins['2mC'],
            'A2': self._origins['2mA'],
            'I2': self._origins['2mI'],
            # ORTHORHOMBIC
            'P 2 2 2': self._origins['222oP'],
            'P 21 2 2': self._origins['222oP'],
            'P 2 21 2': self._origins['222oP'],
            'P 2 2 21': self._origins['222oP'],
            'P 2 21 21': self._origins['222oP'],
            'P 21 2 21': self._origins['222oP'],
            'P 21 21 2': self._origins['222oP'],
            'P 21 21 21': self._origins['222oP'],
            'C 2 2 21': self._origins['222oC'],
            'C 2 2 2': self._origins['222oC'],
            'F 2 2 2': self._origins['222oF'],
            'I 2 2 2': self._origins['222oI'],
            'I 21 21 21': self._origins['222oI'],
            # TETRAGONAL
            'P 4': self._origins['4tP'],
            'P 41': self._origins['4tP'],
            'P 42': self._origins['4tP'],
            'P 43': self._origins['4tP'],
            'I 4': self._origins['4tI'],
            'I 41': self._origins['4tI'],
            'P 4 2 2': self._origins['422tP'],
            'P 4 21 2': self._origins['422tP'],
            'P 41 2 2': self._origins['422tP'],
            'P 41 21 2': self._origins['422tP'],
            'P 42 2 2': self._origins['422tP'],
            'P 42 21 2': self._origins['422tP'],
            'P 43 2 2': self._origins['422tP'],
            'P 43 21 2': self._origins['422tP'],
            'I 4 2 2': self._origins['422tI'],
            'I 41 2 2': self._origins['422tI'],
            # TRIGONAL
            'P 3': self._origins['3hP'],
            'P 31': self._origins['3hP'],
            'P 32': self._origins['3hP'],
            'H 3': self._origins['3hR_1'],
            'R 3': self._origins['3hR_2'],
            'P 3 1 2': self._origins['312hP'],
            'P 31 1 2': self._origins['312hP'],
            'P 32 1 2': self._origins['312hP'],
            'P 3 2 1': self._origins['321hP'],
            'P 31 2 1': self._origins['321hP'],
            'P 32 2 1': self._origins['321hP'],
            'H 3 2': self._origins['32hR_1'],
            'R 3 2': self._origins['32hR_2'],
            # HEXAGONAL
            'P 6': self._origins['6hP'],
            'P 61': self._origins['6hP'],
            'P 65': self._origins['6hP'],
            'P 62': self._origins['6hP'],
            'P 64': self._origins['6hP'],
            'P 63': self._origins['6hP'],
            'P 6 2 2': self._origins['622hP'],
            'P 61 2 2': self._origins['622hP'],
            'P 65 2 2': self._origins['622hP'],
            'P 62 2 2': self._origins['622hP'],
            'P 64 2 2': self._origins['622hP'],
            'P 63 2 2': self._origins['622hP'],
            # CUBIC
            'P 2 3': self._origins['23cP'],
            'P 21 3': self._origins['23cP'],
            'F 2 3': self._origins['23cF'],
            'I 2 3': self._origins['23cI'],
            'I 21 3': self._origins['23cI'],
            'P 4 3 2': self._origins['432cP'],
            'P 42 3 2': self._origins['432cP'],
            'P 43 3 2': self._origins['432cP'],
            'P 41 3 2': self._origins['432cP'],
            'F 4 3 2': self._origins['432cF'],
            'F 41 3 2': self._origins['432cF'],
            'I 4 3 2': self._origins['432cI'],
            'I 41 3 2': self._origins['432cI'],
        }

        return

    def spaceGroup(self):
        return self._spaceGroup

    def isFloating(self, spaceGroupLabel=None):
        if spaceGroupLabel is not None and self.spaceGroup() != spaceGroupLabel:
            self._getAlternateOrigins(spaceGroupLabel)
        return self._floating

    def redundantAlternateOrigins(self, spaceGroupLabel=None):
        if spaceGroupLabel is not None and self.spaceGroup() != spaceGroupLabel:
            self._getAlternateOrigins(spaceGroupLabel)
        return copy.copy(self._redundantSet)

    def nonRedundantAlternateOrigins(self, spaceGroupLabel=None):
        if spaceGroupLabel is not None and self.spaceGroup() != spaceGroupLabel:
            self._getAlternateOrigins(spaceGroupLabel)
        return copy.copy(self._nonRedundantSet)

    def _getAlternateOrigins(self, spaceGroupLabel):
        """Given a space group label, return a list of (non-redundant) alternate
        origins as a list of float triples"""

        label = spaceGroupLabel
        if label not in self._spacegroup2origin:
            label = self._altlabel(label)

        self._spaceGroup = label
        originl = self._spacegroup2origin[label]

        # We build up a list of the full set (redundant) and also the non-redundant that are
        # the only ones we need to loop through when we are checking
        self._nonRedundantSet = []
        self._redundantSet = []
        self._floating = False
        for o in originl:
            if o[1]:
                self._nonRedundantSet.append(o[0])
            self._redundantSet.append(o[0])

        self._floating = any(map(lambda o: 'x' in o or 'y' in o or 'z' in o, self._redundantSet))
        return

    # symoplib = "/Applications/ccp4-6.4.0/lib/data/symop.lib"
    def _altlabel(self, spaceGroup, symoplib=None):

        if not symoplib:
            symoplib = os.path.join(os.environ['CCP4'], "lib/data/symop.lib")

        for line in open(symoplib, 'r'):
            if "'" in line:
                # Assume first single-quote enclosed string is the one we want
                i = line.index("'")
                j = line.index("'", i + 1)
                sg = line[i + 1 : j]
                if spaceGroup == sg:
                    return line.split()[3]

        raise KeyError(spaceGroup)


class CrystalInfo(object):
    def __init__(self, line=None):
        """foo"""

        self._reset()

        if line:
            self.fromLine(line)

        return

    def _reset(self):

        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.spaceGroup = None
        self.z = None

        return

    def fromLine(self, line):

        self.a = float(line[6:15].strip())
        self.b = float(line[15:24].strip())
        self.c = float(line[24:33].strip())
        self.alpha = float(line[33:40])
        self.beta = float(line[40:47])
        self.gamma = float(line[47:54])
        self.spaceGroup = line[55:66].strip()
        try:
            self.z = int(line[66:70])
        except ValueError:
            # Z-info could be missing (shelxe output pdb)
            pass

        return


class PdbInfo(object):
    """A class to hold information extracted from a PDB file"""

    def __init__(self):

        self.models = []  # List of PdbModel objects

        self.pdbCode = None
        self.title = None  # First line of the title
        self.resolution = None

        # http://www.wwpdb.org/documentation/format33/remarks1.html#REMARK%20280
        self.solventContent = None
        self.matthewsCoefficient = None

        self.crystalInfo = None

        return

    def getSequence(self):
        """Return the sequence for the first model/chain"""

        assert len(self.models) >= 1, "Need at least one model!"
        assert len(self.models[0].chains) >= 1, "Need at least one chain!"
        return self.sequences[0]

    def numAtoms(self, modelIdx=0):
        """Return the total number of ATOM atoms in the model"""
        assert len(self.models) >= 1, "Need at least one model!"
        assert len(self.models[modelIdx].chains) >= 1, "Need at least one chain!"

        natoms = 0
        for chainAtoms in self.models[modelIdx].atoms:
            natoms += len(chainAtoms)

        return natoms

    def numChains(self, modelIdx=0):
        """Return the total number of chains in the model"""
        assert len(self.models) >= 1, "Need at least one model!"
        assert len(self.models[modelIdx].chains) >= 1, "Need at least one chain!"
        return len(self.models[modelIdx].chains)

    def numCalpha(self, modelIdx=0):
        """Return the total number of CA ATOM atoms in the model"""
        assert len(self.models) >= 1, "Need at least one model!"
        assert len(self.models[modelIdx].chains) >= 1, "Need at least one chain!"

        ncalpha = 0
        for chainAtoms in self.models[modelIdx].atoms:
            for atom in chainAtoms:
                if atom.name.strip() == 'CA':
                    ncalpha += 1

        return ncalpha


class PdbModel(object):
    """A class to hold information on a single model in a PDB file"""

    def __init__(self):

        self.pdb = None
        self.serial = None
        self.chains = []  # Ordered list of chain IDs
        self.atoms = []  # List of atoms in each chain

        self.resSeqs = []  # Ordered list of list of resSeqs for each chain - matches order in self.chains
        self.sequences = []  # Ordered list of list of sequences for each chain - matches order in self.chains
        self.caMask = []  # Ordered list of list of booleans of residues with no CA atoms - matches order in self.chains
        self.bbMask = (
            []
        )  # Ordered list of list of boleans of residues with no backbone atoms - matches order in self.chains

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
            self.fromLine(line)

        return

    def _setAtomType(self):
        """This gets overridden in HETATM - otherwise everything the same"""
        self._atomType = "ATOM  "
        return

    def _reset(self):

        self.line = None  # the line we were created from
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

    def _readCharge(self, line):
        s = line[78:80]
        minus = '-'
        signs = ['+', minus]
        mult = +1
        if s[0] in signs:
            sign = s[0]
            val = s[1]
        elif s[1] in signs:
            sign = s[1]
            val = s[0]
        else:
            raise RuntimeError("Error getting charge sign ({0}) from line: {1}".format(line[78:80], line))
        if sign == minus:
            mult = -1
        try:
            return int(val) * mult
        except:
            raise RuntimeError("Error getting charge ({0}) from line: {1}".format(line[78:80], line))

    def _sanityCheck(self, line):
        assert line[0:6] == self._atomType, "Line did not begin with an {0} record!: {1}".format(self._atomType, line)
        assert len(line) >= 54, "Line length was: {0}\n{1}".format(len(line), line)

    def fromLine(self, line):
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
            self.charge = self._readCharge(line)
        return

    def toLine(self):
        """Create a line suitable for printing to a PDB file"""

        s = self._atomType  # 1-6
        s += "{0:5d}".format(self.serial)  # 7-11
        s += " "  # 12 blank
        if len(self.name) != 4:
            raise RuntimeError("Name must be 4 characters long!")
        s += "{0:4}".format(self.name)  # 13-16
        if not self.altLoc:  # 17
            s += " "
        else:
            s += "{0:1}".format(self.altLoc)
        s += "{0:3}".format(self.resName)  # 18-20
        s += " "  # 21 blank
        if not self.chainID:  # 22
            s += " "
        else:
            s += "{0:1}".format(self.chainID)
        s += "{0:4}".format(self.resSeq)  # 23-26
        if not self.iCode:  # 27
            s += " "
        else:
            s += "{0:1}".format(self.iCode)
        s += "   "  # 28-30 blank
        s += "{0:8.3F}".format(self.x)  # 31-38
        s += "{0:8.3F}".format(self.y)  # 39-46
        s += "{0:8.3F}".format(self.z)  # 47-54
        if not self.occupancy:  # 55-60
            s += "      "
        else:
            s += "{0:6.2F}".format(self.occupancy)
        if not self.tempFactor:  # 61-66
            s += "      "
        else:
            s += "{0:6.2F}".format(self.tempFactor)
        s += "      "  # 67-72 blank
        if not self.segID:  # 73-76
            s += "    "
        else:
            s += "{0:>4}".format(self.segID)
        if not self.element:  # 77-78
            s += "  "
        else:
            s += "{0:>2}".format(self.element)
        if not self.charge:  # 79-80
            s += "  "
        else:
            s += "{0:2d}".format(self.charge)
        return s

    def fromHetatm(self, hetatm):
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
            if not slot.startswith("__") and not (
                isinstance(attr, types.MethodType) or isinstance(attr, types.FunctionType)
            ):
                me[slot] = attr
        return "{0} : {1}".format(self.__repr__(), str(me))


class PdbHetatm(PdbAtom):
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

        self.fromLine(line)

    def _reset(self):

        self.idCode = None
        self.resName = None
        self.chainID = None
        self.seqNum = None
        self.iCode = None
        self.stdRes = None
        self.comment = None

        return

    def fromLine(self, line):
        """Initialise from the line from a PDB"""

        assert line[0:6] == "MODRES", "Line did not begin with an MODRES record!: {0}".format(line)

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

        s = "MODRES"  # 1-6
        s += " "  # 7 blank
        s += "{0:4}".format(self.idCode)  # 8-11
        s += " "  # 12 blank
        s += "{0:>3}".format(self.resName)  # 13-15
        s += " "  # 16 blank
        if not self.chainID:  # 17
            s += " "
        else:
            s += "{0:1}".format(self.chainID)
        s += " "  # 18 blank
        s += "{0:4d}".format(self.seqNum)  # 19-22
        if not self.iCode:  # 23
            s += " "
        else:
            s += "{0:1}".format(self.iCode)
        s += " "  # 24 blank
        s += "{0:>3}".format(self.stdRes)  # 25-27
        s += "  "  # 28-29 blank
        if self.comment:  # 30-70
            s += "{:<}".format(self.comment)

        return s

    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not (
                isinstance(attr, types.MethodType) or isinstance(attr, types.FunctionType)
            ):
                me[slot] = attr

        return "{0} : {1}".format(self.__repr__(), str(me))
