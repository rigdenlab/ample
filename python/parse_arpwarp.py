'''
Created on 28 Nov 2013

@author: jmht
'''

import unittest

class ArpwarpLogParser(object):
    """
    Class to mine information from an arpwarp log
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

        #CAPTURE=False
        cycle=0
        res_built=0
        Rfact=1.0
        Rfree=1.0
        with open( logfile, 'r') as fh:
            line = fh.readline()
            while line:
                if "Building Cycle" in line and "Atomic shape" in line:
                    cycle=cycle+1
                if "Estimated correctness" in line:
                    pass
                if "docked in sequence" in line:
                    res_built=int(line.split()[2].replace("(",""))
                if "After refmac" in line:
                    #CAPTURE = True
                    Rfact=float( line.split("R =")[1].split()[0] )
                    Rfree=float(line.split("Rfree =")[1].split()[0].replace(").","").replace(")",""))
                    #Rcycle=line.split()[1]
                    if cycle==0:
                        # jmht - not sure this correct - shouldn't we be using
                        #  Starting model:  R = 0.4486 (Rfree = 0.480). 
                        self.initRfact=Rfact
                        self.initRfree=Rfree
                #if CAPTURE and "------------------" in line:
                #    CAPTURE = False
    
                line = fh.readline()
            #End while
    
        self.res_built=res_built
        self.finalRfact=Rfact
        self.finalRfree=Rfree
        
        return
    
class Test(unittest.TestCase):

    def testParse1(self):
        """parse 2bhw"""
        
        logFile = "/opt/ample-dev1/tests/testfiles/arpwarp.log"
        ap = ArpwarpLogParser( logFile )
        
        #self.assertEqual( ap.initRfact, 0.452)
        #self.assertEqual( ap.initRfree, 0.4199)
        self.assertEqual( ap.finalRfree, 0.242)
        self.assertEqual( ap.finalRfact, 0.2157)
        self.assertEqual( ap.res_built, 44)
        
        return


def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testParse1'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())

