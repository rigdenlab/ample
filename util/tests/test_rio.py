
import os
import unittest
from ample.util import rio
from ample.testing import constants

class TestContacts( unittest.TestCase ):
    
    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_dir = constants.AMPLE_DIR
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')

    def test_parse1(self):
        logfile = os.path.join( self.testfiles_dir, "ncont1.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 26 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 0 )
        self.assertEqual( contactData.rioBackwards, 0 )
    
    def test_parse2(self):
        logfile = os.path.join( self.testfiles_dir, "ncont2.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 10 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 7 )
        self.assertEqual( contactData.rioBackwards, 7 )
    
    def test_parse3(self):
        logfile = os.path.join( self.testfiles_dir, "ncont3.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 14 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 10 )
        self.assertEqual( contactData.rioBackwards, 0 )
    
    def test_parse4(self):
        logfile = os.path.join( self.testfiles_dir, "ncont4.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 56 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 55 )
        self.assertEqual( contactData.rioBackwards, 0 )
    
    def test_parse5(self):
        logfile = os.path.join( self.testfiles_dir, "ncont5.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 77 )
        self.assertEqual( contactData.rioInRegister, 19 )
        self.assertEqual( contactData.rioOoRegister, 54 )
        self.assertEqual( contactData.rioBackwards,16 )

    def test_parse7(self):
        logfile = os.path.join( self.testfiles_dir, "ncont7.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 18 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 0 )
        self.assertEqual( contactData.rioBackwards,0 )
    
    def test_parse8(self):
        logfile = os.path.join( self.testfiles_dir, "ncont8.log" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )
        c.analyseRio(contactData)
        self.assertEqual( contactData.numContacts, 9 )
        self.assertEqual( contactData.rioInRegister, 0 )
        self.assertEqual( contactData.rioOoRegister, 0 )
        self.assertEqual( contactData.rioBackwards,0 )


    def test_helix5(self):
        logfile = os.path.join( self.testfiles_dir, "ncont5.log" )
        dssplog = os.path.join( self.testfiles_dir, "3RA3.dssp" )
        c = rio.Rio()
        contactData = rio.RioData()
        c.parseNcontLog( contactData, logfile=logfile, clean_up=False )       
        sequence = c.helixFromContacts( contactData.contacts, dssplog )
        self.assertEqual( "NARLKQEIAALEYEIAAL", sequence )

if __name__ == "__main__":
    unittest.main()
