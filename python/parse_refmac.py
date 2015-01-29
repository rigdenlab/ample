'''
Created on 28 Nov 2013

@author: jmht
'''

import os
import unittest

class RefmacLogParser(object):
    """
    Class to mine information from a refmac log
    """
    
    def __init__(self,logfile=None):
        
        if logfile:
            self.parse( logfile )
        return
        
    def parse(self, logfile ):
        """Parse 
        """
        #self.initRfree=1.0
        #self.initRfact=1.0
        self.finalRfree=1.0
        self.finalRfact=1.0

        Rfact=1.0
        Rfree=1.0
        with open( logfile, 'r') as fh:
            line = fh.readline()
            # Brain-dead approach...
            while line:
                line=line.strip()
                if line.startswith("R factor"):
                    Rfact=float(line.split()[3])
                if line.startswith("R free"):
                    Rfree=float(line.split()[3])
    
                line = fh.readline()
            #End while
    
        self.finalRfact=Rfact
        self.finalRfree=Rfree
        return
    
class Test(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return
    
    def testParse1(self):
        """parse 2bhw"""
        
        logFile = os.path.join(self.testfiles_dir,"refmac.log")
        ap = RefmacLogParser( logFile )
        
        #self.assertEqual( ap.initRfact, 0.452)
        #self.assertEqual( ap.initRfree, 0.4199)
        self.assertEqual( ap.finalRfact, 0.5207)
        self.assertEqual( ap.finalRfree, 0.5289)
        
        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testParse1'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())

