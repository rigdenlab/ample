#!/usr/bin/env python

import ample_util

class Csymmatch( object ):
    
    def __init__(self):
        
        self.logfile=None
        
        return
  
    def run( self, refPdb=None, inPdb=None, outPdb=None, connectivityRadius=None, originHand=True ):
        """FOO
        """
        
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
    
    def origin( self,  logfile=None ):
        """Return the change of origin"""
        
        if not logfile:
            logfile = self.logfile
        
        for line in open( logfile, 'r' ):
            if "Change of hand" in line:
                # For our work, if there is a change of hand it didn't work
                fields = line.split()
                if fields[ 4 ] == "YES":
                    return False
            
            if "Change of origin:" in line:
                # Find position of brackets
                i1 = line.index( "(" )
                i2 = line.index( ")" )
                oline = line[ i1+1:i2 ]
                x,y,z = oline.split(",")
                return [ float(x), float(y), float(z) ]
        
        return False

if __name__ == "__main__":
    logfile = "/home/jmht/Documents/test/ncont/new/1GMJ/COH/shelxe_phaser_loc0_ALL_poly_ala_trunc_1.732018_rad_2_UNMOD_csymmatch.pdb.log"
    c = Csymmatch()
    print c.origin( logfile=logfile )


