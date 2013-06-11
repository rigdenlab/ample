'''
Created on 28 May 2013

@author: jmht
'''

import ample_util
import logging
import os


class CifParser(object):
    """Class for manipulating CIF files."""
    
    
    def sfcif2mtz(self, cifpath, mtzpath=None ):
        """Convert a CIF containing structure factors to an MTZ file."""
        
        if not mtzpath:
            name = os.path.basename(cifpath)
            mtzpath = os.path.join( os.getcwd(), name+".mtz" )
            
        logfile = "cif2mtz.log"
        
        cmd = [ "cif2mtz", "hklin", cifpath, "hklout", mtzpath ]
        
        # Empty stdin - just send EOF
        stdin=""
        
        retcode = ample_util.run_command(cmd, logfile=logfile, stdin=stdin)
        
        if retcode != 0:
            msg = "Error running sfcif2mtz. Check the logfile: {0}".format(logfile)
            logging.critical(msg)
            raise RuntimeError, msg
        
        return mtzpath
        


if __name__ == '__main__':
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    cifpath="/media/data2/Shared/EVFOLD/3B45/3b45-sf.cif"
    
    cp = CifParser()
    cp.sfcif2mtz(cifpath, mtzpath="foo.mtz")
