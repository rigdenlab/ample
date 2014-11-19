'''
Created on 28 Nov 2013

@author: jmht
'''
import os
import unittest

class BuccaneerLogParser(object):
    """
    Class to mine information from a buccaneer log
    """
    
    def __init__(self,logfile=None):
        
        if logfile:
            self.parse( logfile )
        return
        
    def parse(self, logfile ):
        """Parse 
        """
    
        self.initRfree=1.0
        self.finalRfree=1.0
        self.initRfact=1.0
        self.finalRfact=1.0
    
        fh = open( logfile, 'r')
        line = fh.readline()
        while line:
            if line.startswith("           R factor"):
                self.initRfact = float( line.split()[2] )
                self.finalRfact = float( line.split()[3] )
            if "R free" in line:
                self.initRfree = float( line.split()[2] )
                self.finalRfree = float( line.split()[3] )

            line = fh.readline()
        #End while
        
        fh.close()
        
        return
    
class Test(unittest.TestCase):


    def setUp(self):
        
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ampleDir = os.sep.join( paths[ : -1 ] )
        self.testfilesDir = os.sep.join( paths[ : -1 ] + [ 'tests', 'testfiles' ] )
        
        return

    def testParse1(self):
        """parse 2bhw"""
        
        logFile = os.path.join(self.testfilesDir,"buccaneer.log")
        bp = BuccaneerLogParser( logFile )
        
        self.assertEqual( bp.initRfree, 0.3781)
        self.assertEqual( bp.finalRfree, 0.3222)
        self.assertEqual( bp.initRfact, 0.3708)
        self.assertEqual( bp.finalRfact, 0.2979)
        
        return
    
def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testParse1'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
