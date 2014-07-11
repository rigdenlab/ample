'''
Created on 28 Nov 2013

@author: jmht
'''

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

    def testParse1(self):
        """parse 2bhw"""
        
        logFile = "/Users/jmht/Documents/AMPLE/data/buccaneer.log"
        bp = BuccaneerLogParser( logFile )
        
        self.assertEqual( bp.initRfree, 0.5642)
        self.assertEqual( bp.finalRfree, 0.2841)
        self.assertEqual( bp.initRfact, 0.5548)
        self.assertEqual( bp.finalRfact, 0.2659)
        
        return
    
if __name__ == "__main__":
    
    unittest.main()