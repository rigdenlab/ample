#!/usr/bin/python2.6

#

import glob
import os
import logging

import ample_util

class Scwrl( object ):
    
    def __init__(self, scwrl_exe=None, workdir=None ):
        self.workdir = workdir
        if self.workdir is None: self.workdir = os.getcwd()
        if not ample_util.is_exe(scwrl_exe): raise RuntimeError("scwrl_exe {0} cannot be found.".format(scwrl_exe))
        self.scwrl_exe = scwrl_exe
        return
    
    def add_sidechains(self, pdbin=None, pdbout=None, sequence=None, hydrogens=False):
        """Add the specified sidechains to the pdb"""
        
        cmd = [ self.scwrl_exe, "-i", pdbin, "-o", pdbout ]
        
        # Not needed by default
        if sequence is not None:
            sequenceFile = os.path.join( self.workdir, "sequence.file")
            with open( sequenceFile, 'w' ) as w:
                w.write( sequence + os.linesep )
            cmd += [ "-s",  sequenceFile ]
        
        # Don't output hydrogens
        if not hydrogens: cmd += ['-h']
            
        logfile = os.path.abspath("scwrl.log")
        retcode = ample_util.run_command(cmd, logfile=logfile)
        
        if retcode != 0:
            raise RuntimeError,"Error running Scwrl - please check the logfile: {0}".format(logfile)
        else:
            os.unlink(logfile)
        
        return
    
    def process_directory(self, in_dir, out_dir, prefix="scwrl" ):
        logging.info('Adding sidechains with SCWRL to models in directory: {0}'.format(in_dir))
        count = 0
        for pdb in glob.glob( os.path.join( in_dir, '*.pdb') ):
            pdbout = ample_util.filename_append( pdb, prefix, directory=out_dir )
            self.add_sidechains(pdbin=pdb, pdbout=pdbout)
            count += 1
        
        logging.info('Processed {0} models with SCWRL into directory: {1}'.format(count, out_dir))
        return 
