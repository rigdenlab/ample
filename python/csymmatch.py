#!/usr/bin/env python

import os
import unittest

# Our imports
import ample_util
import pdb_edit

class Csymmatch( object ):
    
    def __init__(self):
        
        self._reset()
        return
    

    def _reset(self):
        self.logfile=None
        self.parsed=False
        self.changeOfHand = False
        self.changeOfOrigin = None
        self.chainShifts = {}
        return
        
    def run( self, refPdb=None, inPdb=None, outPdb=None, connectivityRadius=None, originHand=True ):
        """FOO
        """
        
        self._reset()
        
        self.logfile = outPdb +".log"
        cmd= [ 'csymmatch',
              "-pdbin-ref",
              refPdb,
              "-pdbin",
              inPdb,
              "-pdbout",
              outPdb
              ]
        
        if originHand:
            cmd += [  "-origin-hand" ]

        if connectivityRadius:
            cmd += [ "-connectivity-radius", connectivityRadius ]

        retcode = ample_util.run_command(cmd=cmd, logfile=self.logfile, dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        return
    
    def parseLog( self,  logfile=None ):
        """Parse the log"""
        
        if logfile is None:
            logfile = self.logfile
        
        assert logfile
        
        capturing=0
        currentChain=None
        self.chainShifts = {}
        for line in open( logfile, 'r' ):
            
            line = line.strip()
            if not line:
                continue
            
            if line.startswith("Chain:"):
                capturing=1
                f = line.split()
                currentChain=f[1]
                resStart = int( f[2][:-1] )
                resEnd = int( f[3] )
                if currentChain not in self.chainShifts:
                    self.chainShifts[ currentChain ] = []
                self.chainShifts[ currentChain ].append( { 'resStart' : resStart, 'resEnd' : resEnd } )
                continue
            
            if capturing == 1 or capturing == 2:
                capturing += 1 # Symmetry Opterator & Lattice Shift
            elif capturing == 3:
                capturing += 1
                score = float( line.split()[3] )
                self.chainShifts[ currentChain ][ -1 ][ 'score' ] = score
            elif capturing == 4:
                # Singnifies end of chain blocks
                continue
                
            if "Change of hand" in line:
                # For our work, if there is a change of hand it didn't work
                fields = line.split()
                if fields[ 4 ] == "YES":
                    self.changeOfHand = True
            
            if "Change of origin:" in line:
                # Find position of brackets
                i1 = line.index( "(" )
                i2 = line.index( ")" )
                oline = line[ i1+1:i2 ]
                x,y,z = oline.split(",")
                self.changeOfOrigin = [ float(x), float(y), float(z) ]
        
        self.parsed=True
        return

    
    def origin( self,  logfile=None, failOnChangeOfHand=True ):
        """Return the change of origin.
        Csymmatch will always return something so we use a changeOfHand as indication of failure
        """
        if not self.parsed:
            self.parseLog()
        if failOnChangeOfHand and self.changeOfHand:
            return False
        return self.changeOfOrigin
    
    def averageScore(self):
        if not self.parsed:
            self.parseLog()
            
        if not len( self.chainShifts ):
            return False
        
        score=0.0
        count=0
        for chain in self.chainShifts.keys():
            for shift in self.chainShifts[ chain ]:
                score += shift['score']
                count += 1
        
        return score/count
    
    def wrapModelToNative(self, mrPdb, nativePdb, origin=[0.0,0.0,0.0], csymmatchPdb=None, workdir=None ):
        """Take a pdb and wrap it onto the nativePdb using csymmatch.
        If origin is not [0.0,0.0,0.0] we also move the structure onto the new origin before wrapping"""
        
        if workdir is None:
            workdir = os.getcwd()
        
        assert os.path.isfile( mrPdb ) and os.path.isfile( nativePdb )
        
        if origin != [ 0.0, 0.0, 0.0 ]:
            # Move pdb to new origin
            #ostr="origin{0}".format(i)
            ostr="o{0}".format( origin ).replace(" ","" )
            originMrPdb = ample_util.filename_append( filename=mrPdb, astr=ostr, directory=workdir )
            pdb_edit.PDBEdit().translate( inpdb=mrPdb, outpdb=originMrPdb, ftranslate=origin )
            mrPdb = originMrPdb
        
        if csymmatchPdb is None:
            csymmatchPdb = ample_util.filename_append( filename=mrPdb,
                                                       astr="csymmatch",
                                                       directory=workdir )
        self.run( refPdb=nativePdb,
                  inPdb=mrPdb,
                  outPdb=csymmatchPdb,
                  originHand=False )
        
        if not os.path.isfile( csymmatchPdb ):
            raise RuntimeError,"Error generating csymmatchPdb" 
               
        return csymmatchPdb
        

class TestContacts( unittest.TestCase ):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def testParse1(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        logfile = os.path.join( self.testfiles_dir, "csymmatch1.log" )
        
        c = Csymmatch()
        c.parseLog( logfile=logfile )

        self.assertEqual( c.changeOfHand, False )
        self.assertEqual( c.changeOfOrigin, None )
        self.assertEqual( c.chainShifts, {'a': [{'resStart': 14, 'score': 0.403697, 'resEnd': 14}, {'resStart': 17, 'score': 0.247688, 'resEnd': 17}, {'resStart': 32, 'score': 0.528113, 'resEnd': 44}, {'resStart': 49, 'score': 0.268943, 'resEnd': 51}]}
 )
        self.assertEqual( c.averageScore(), 0.36211025)
        return
    
    def testParse2(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        
        logfile = os.path.join( self.testfiles_dir, "csymmatch2.log" )
        
        c = Csymmatch()
        c.parseLog( logfile=logfile )

        self.assertEqual( c.changeOfHand, True )
        self.assertEqual( c.changeOfOrigin, [ 0, 0.5625, 0 ] )
        self.assertEqual( c.chainShifts, {'A': [{'resStart': 1, 'score': 0.440815, 'resEnd': 27} ],'B': [{'resStart': 1, 'score': 0.558538, 'resEnd': 6} ] }
 )
        self.assertEqual( c.averageScore(), 0.49967649999999997)
        return
    

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(TestContacts('testParse1'))
    suite.addTest(TestContacts('testParse2'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
