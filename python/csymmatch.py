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
        
    def run(self, refPdb=None, inPdb=None, outPdb=None, connectivityRadius=None, originHand=True, cleanup=False):
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
        
        if retcode != 0: raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        if cleanup: os.unlink(self.logfile)
        return
    
    def parseLog(self,  logfile=None, cleanup=True):
        """Parse the log"""
        
        if logfile is None: logfile = self.logfile
        assert logfile
        
        capturing=0
        currentChain=None
        self.chainShifts = {}
        for line in open( logfile, 'r' ):
            
            line = line.strip()
            if not line: continue
            
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
        if cleanup: os.unlink(logfile)
        return
    
    def origin(self, logfile=None, failOnChangeOfHand=True):
        """Return the change of origin.
        Csymmatch will always return something so we use a changeOfHand as indication of failure
        """
        if not self.parsed: self.parseLog()
        if failOnChangeOfHand and self.changeOfHand: return False
        return self.changeOfOrigin
    
    def averageScore(self):
        if not self.parsed: self.parseLog()
            
        if not len(self.chainShifts): return False
        
        score=0.0
        count=0
        for chain in self.chainShifts.keys():
            for shift in self.chainShifts[ chain ]:
                score += shift['score']
                count += 1
        
        return score/count
    
    def wrapModelToNative(self,
                          mrPdb,
                          nativePdb,
                          origin=[0.0,0.0,0.0],
                          csymmatchPdb=None,
                          workdir=None,
                          cleanup=True):
        """Take a pdb and wrap it onto the nativePdb using csymmatch.
        If origin is not [0.0,0.0,0.0] we also move the structure onto the new origin before wrapping"""
        
        if workdir is None: workdir = os.getcwd()
        
        assert os.path.isfile(mrPdb) and os.path.isfile(nativePdb),"Cannot find: {0} or {1}".format(mrPdb,nativePdb)
        
        originMrPdb = None
        if origin != [ 0.0, 0.0, 0.0 ]:
            # Move pdb to new origin
            #ostr="origin{0}".format(i)
            ostr="o{0}".format( origin ).replace(" ","" )
            originMrPdb = ample_util.filename_append(filename=mrPdb, astr=ostr, directory=workdir )
            pdb_edit.translate(inpdb=mrPdb, outpdb=originMrPdb, ftranslate=origin)
            mrPdb = originMrPdb
        
        if csymmatchPdb is None:
            csymmatchPdb = ample_util.filename_append(filename=mrPdb, astr="csymmatch", directory=workdir)
        
        self.run(refPdb=nativePdb,
                 inPdb=mrPdb,
                 outPdb=csymmatchPdb,
                 originHand=False,
                 cleanup=cleanup)
        
        if not os.path.isfile( csymmatchPdb ): raise RuntimeError,"Error generating csymmatchPdb"
        
        if cleanup and originMrPdb: os.unlink(originMrPdb)
               
        return csymmatchPdb

