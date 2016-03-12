
import unittest
from ample.util import pdb_model

class Test(unittest.TestCase):

    def test_read_atom(self):
        #See if we can read an atom line
        line = "ATOM     41  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N  "
        a = pdb_model.PdbAtom( line )
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
    
    def test_read_atom2(self):
        #Round-trip an atom line
        line = "ATOM     28  C   ALA A  12     -27.804  -2.987  10.849  1.00 11.75      AA-- C "
        a = pdb_model.PdbAtom( line )
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
    
    def test_read_atom3(self):
        #Round-trip an atom line
        line = "ATOM    160  NH1 ARG A  21      57.124  31.377  40.357  1.00 35.50           N1+"
        a = pdb_model.PdbAtom( line )
        self.assertEqual(a.serial,160)
        self.assertEqual(a.name,' NH1')
        self.assertEqual(a.altLoc,None)
        self.assertEqual(a.resName,'ARG')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,21)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,57.124)
        self.assertEqual(a.y,31.377)
        self.assertEqual(a.z,40.357)
        self.assertEqual(a.occupancy,1.00)
        self.assertEqual(a.tempFactor,35.50)       
        self.assertEqual(a.element,'N')
        self.assertEqual(a.charge,1)
    
    def test_read_atom4(self):
        #Round-trip an atom line
        line = "ATOM    183  OD2 ASP A  24      70.534  30.495  41.026  1.00 35.00           O1-"
        a = pdb_model.PdbAtom( line )
        self.assertEqual(a.serial,183)
        self.assertEqual(a.name,' OD2')
        self.assertEqual(a.altLoc,None)
        self.assertEqual(a.resName,'ASP')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.resSeq,24)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.x,70.534)
        self.assertEqual(a.y,30.495)
        self.assertEqual(a.z,41.026)
        self.assertEqual(a.occupancy,1.00)
        self.assertEqual(a.tempFactor,35.00)       
        self.assertEqual(a.element,'O')
        self.assertEqual(a.charge,-1)
           
    def test_write_atom1(self):
        #Round-trip an atom line
        line = "ATOM     41  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N  "
        a = pdb_model.PdbAtom( line )
        self.assertEqual( a.toLine(), line )
           
    def test_read_hetatm(self):
        #See if we can read a hetatom line
        line = "HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG  "
        a = pdb_model.PdbHetatm( line )
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
    
    def test_write_hetatm(self):
        #Round-trip an atom line        
        line = "HETATM 8239   O1 SO4 A2001      11.191 -14.833 -15.531  1.00 50.12           O  "
        a = pdb_model.PdbHetatm( line )
        self.assertEqual( a.toLine(), line )
  
    def test_write_hetatm2(self):
        #Round-trip an atom line
        line = "HETATM    7  SD  FME A   1     -60.099  -1.874   3.446  1.00216.81           S  "
        a = pdb_model.PdbHetatm( line )
        self.assertEqual( a.toLine(), line )
   
    def test_read_modres(self):
        #See if we can read a modres line
        line = "MODRES 1IL2 1MG D 1937    G  1N-METHYLGUANOSINE-5'-MONOPHOSPHATE"
        a = pdb_model.PdbModres( line )
        self.assertEqual(a.idCode,"1IL2")
        self.assertEqual(a.resName,'1MG')
        self.assertEqual(a.chainID,'D')
        self.assertEqual(a.seqNum,1937)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.stdRes,'G')
        self.assertEqual(a.comment,"1N-METHYLGUANOSINE-5'-MONOPHOSPHATE")
    
    def test_write_modres(self):
        #Round-trip a modres line
        line = "MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE"
        a = pdb_model.PdbModres( line )
        self.assertEqual(a.idCode,"2R0L")
        self.assertEqual(a.resName,'ASN')
        self.assertEqual(a.chainID,'A')
        self.assertEqual(a.seqNum,74)
        self.assertEqual(a.iCode,None)
        self.assertEqual(a.stdRes,'ASN')
        self.assertEqual(a.comment,"GLYCOSYLATION SITE")
        self.assertEqual( a.toLine(), line )
    
    def test_read_crystal_info(self):
        # See if we can read a cryst1 line
        line = "CRYST1  117.000   15.000   39.000  90.00  90.00  90.00 P 21 21 21    8 "
        a = pdb_model.CrystalInfo( line )
        self.assertEqual(a.a,117.000)
        self.assertEqual(a.b,15.000)
        self.assertEqual(a.c,39.000)
        self.assertEqual(a.alpha,90)
        self.assertEqual(a.beta,90)
        self.assertEqual(a.gamma,90)
        self.assertEqual(a.spaceGroup,"P 21 21 21")
        self.assertEqual(a.z,8)

if __name__ == "__main__":
    unittest.main()