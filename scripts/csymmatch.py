#!/usr/bin/env python

import os
import sys

import ample_util

class Csymmatch( object ):
    
    def __init__(self):
        pass
  
    def run( self, refPdb=None, inPdb=None, outPdb=None, connectivityRadius=None, originHand=True ):
        """FOO
        """
        
        logfile = outPdb +".log"
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

        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        return
    
    def origin( self,  logfile=None ):
        """Return the change of origin"""
        
        for line in open( logfile, 'r' ):
            if "Change of origin:" in line:
                # Find position of brackets
                i1 = line.index( "(" )
                i2 = line.index( "()" )
                oline = line[ i1:i2 ]
                x,y,z = oline.split(",")
                return [ float(x), float(y), float(z) ]
        
        return False
        
        