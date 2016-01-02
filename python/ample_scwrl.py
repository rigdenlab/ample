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
        return os.path.abspath(pdbout)
    
    def process_directory(self, in_dir, out_dir, prefix="scwrl" ):
        logging.info('Adding sidechains with SCWRL to models in directory: {0}'.format(in_dir))
        self.process_models(glob.glob( os.path.join( in_dir, '*.pdb') ), out_dir, prefix=prefix)
        return
    
    def process_models(self, models, out_dir, prefix="scwrl" ):
        logging.info('Adding sidechains with SCWRL to models')
        out_pdbs = []
        for i, pdb in enumerate(models):
            out_pdbs.append(self.add_sidechains(pdbin=pdb,
                                                pdbout=ample_util.filename_append(pdb,
                                                                                  prefix,
                                                                                  directory=out_dir)))
        logging.info('Processed {0} models with SCWRL into directory: {1}'.format(i+1, out_dir))
        return out_pdbs