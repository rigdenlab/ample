#!/usr/bin/env python

import glob
import os
import re

import ample_util
import pdb_edit
import residue_map 

class MaxclusterData(object):
    
    def __init__(self):
        self.length    = None
        self.pairs     = None
        self.rmsd      = None
        self.maxsub    = None
        self.tm        = None
        self.msi       = None
        self.modelName = None
        self.pdb       = None
        return

    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
            
        return "{0} : {1}".format(self.__repr__(),str(me))

class Maxcluster(object):
    """
        
    # Extract the first chain from the nativePdb
    
    # Create a residueSequenceMap and see if the residues match
    
    # If not use keep_matching to create a nativePdb that has the correct residue sequence`
    
    # Run Maxcluster to compare the models to the native
    """
    
    def __init__(self):

        #self.maxclusterExe = "/Users/jmht/Documents/AMPLE/programs/maxcluster"
        self.maxclusterExe = "/opt/maxcluster/maxcluster"
        
        return
        
    def compareDirectory(self, nativePdbInfo=None, 
                               resSeqMap=None, 
                               modelsDirectory=None,
                               workdir=None ):

        self.data = []
        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
        
        #refModel = os.path.join( modelsDirectory, "S_00000001.pdb" )
        nativePdb = self.prepareNative( nativePdbInfo=nativePdbInfo, resSeqMap=resSeqMap )
        
        logfile = os.path.join( self.workdir, "maxclusterD.log" )
        self.runCompareDirectory( nativePdb=nativePdb, modelsDirectory=modelsDirectory, logfile=logfile )
        self.parseLogDirectory( logfile=logfile )

        return
    
    def compareSingle(self, nativePdb=None, modelPdb=None, sequenceIndependant=True, rmsd=False, workdir=None ):

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
        cmd = [ self.maxclusterExe, "-e", nativePdb, "-p", modelPdb ]
        
        if sequenceIndependant:
            cmd.append( "-in" )
        
        if rmsd:
            cmd.append( "-rmsd" )
        
        logfile = ample_util.filename_append( filename=modelPdb, 
                                    astr="maxcluster", 
                                    directory=self.workdir )
        
        if rmsd:
            logfile = os.path.splitext( logfile )[0] + "_rmsd.log"
        else:
            logfile = os.path.splitext( logfile )[0] + ".log"
        self.maxclusterLogfile = logfile
        
        #print "running cmd "," ".join( cmd )
        retcode = ample_util.run_command( cmd, logfile=self.maxclusterLogfile, dolog=False )
        
        if retcode != 0:
            msg = "non-zero return code for maxcluster in runMaxcluster!"
            #logging.critical( msg )
        
        if rmsd:
            data = self.parseLogSingleRmsd()
        else:
            data = self.parseLogSingleTm()
        
        return data
        
    def prepareNative(self, nativePdbInfo=None, resSeqMap=None ):
        """do stuff"""
        
        # Find out how many chains and extract the first if > 1
        if len( nativePdbInfo.models ) > 1:
            raise RuntimeError,"More than 1 model."
        
        # Check if > 1 chain
        chainID=None
        PE = pdb_edit.PDBEdit()
        if len( nativePdbInfo.models[0].chains ) > 1:
            
            chainID=nativePdbInfo.models[0].chains[0]
        
            # Assume native is standardised
            # Extract the chain if > 1
            nativePdbChain = ample_util.filename_append( filename=nativePdbInfo.pdb,
                                                         astr="chain{0}".format(chainID) )
            PE.extract_chain(nativePdbInfo.pdb, nativePdbChain, chainID)
            nativePdb = nativePdbChain
        else:
            nativePdb = nativePdbInfo.pdb

        if not resSeqMap.resSeqMatch():
            
            # We need to create a copy of the native with numbering matching the model
            nativeRenumber = ample_util.filename_append( filename=nativePdb,
                                                         astr="ren".format(chainID) )
            PE.match_resseq( targetPdb=nativePdb, outPdb=nativeRenumber, resMap=resSeqMap )
            nativePdb = nativeRenumber
        
        return nativePdb
        
    def parseLogDirectory(self, logfile=None ):
        
        self.data = []
        assert logfile
        
        #INFO  : 1000. 2XOV_clean_ren.pdb vs. /media/data/shared/TM/2XOV/models/S_00000444.pdb  Pairs=  36, RMSD= 3.065, MaxSub=0.148, TM=0.192, MSI=0.148
        for line in open( logfile, 'r' ):
            
            if re.match( "INFO *: .* vs\. .* Pairs=", line ):
                
                # Remove spaces after = 
                line = re.sub("= +", "=", line )
                # Now remove commas
                line = line.replace(",","")
                
                fields = line.split()
                
                d = MaxclusterData()
                d.pdb = fields[5]
                d.modelName = os.path.splitext( os.path.basename( fields[5] ) )[0]
                
                label, value = fields[6].split( "=" )
                assert label == "Pairs"
                d.pairs = int( value )
                
                label, value = fields[7].split( "=" )
                assert label == "RMSD"
                d.rmsd = float( value )
                
                label, value = fields[8].split( "=" )
                assert label == "MaxSub"
                d.maxsub = float( value )
                
                label, value = fields[9].split( "=" )
                assert label == "TM"
                d.tm = float( value )
                
                label, value = fields[10].split( "=" )
                assert label == "MSI"
                d.msi = float( value )
                
                self.data.append( d )
                
        return

    def parseLogSingleTm(self, logfile=None):
        
        if logfile is None:
            logfile = self.maxclusterLogfile
        
        assert logfile
        
        d = MaxclusterData()
        for line in open( logfile, 'r' ):
            
            line = line.strip()
            #"Iter 1: Pairs=  14, RMSD= 0.155, MAXSUB=0.855. Len=  15. gRMSD= 0.673, TM=0.858
            if re.match( "Iter \d?: ?Pairs=", line ):
                
                # Remove spaces after = 
                line = re.sub("= +", "=", line )
                
                # Remove trailing . after numbers
                line = re.sub("\. +", " ", line )
                
                # Now remove commas
                line = line.replace(",","")
                
                fields = line.split()
                
                label, value = fields[2].split( "=" )
                assert label == "Pairs"
                d.pairs = int( value )
                
                label, value = fields[3].split( "=" )
                assert label == "RMSD"
                d.rmsd = float( value )
                
                label, value = fields[4].split( "=" )
                assert label == "MAXSUB"
                d.maxsub = float( value )
                
                label, value = fields[5].split( "=" )
                assert label == "Len"
                d.length = float( value )
                
                # skip gRMSD
                
                label, value = fields[7].split( "=" )
                assert label == "TM"
                d.tm = float( value )
        
        return d

    def parseLogSingleRmsd(self, logfile=None):
        
        if logfile is None:
            logfile = self.maxclusterLogfile
        
        assert logfile
        
        d = MaxclusterData()
        
        #INFO  : 1000. 2XOV_clean_ren.pdb vs. /media/data/shared/TM/2XOV/models/S_00000444.pdb  Pairs=  36, RMSD= 3.065, MaxSub=0.148, TM=0.192, MSI=0.148
        for line in open( logfile, 'r' ):
            
            line = line.strip()
            
            if line.startswith("RMSD="):
                # RMSD= 0.132 (Pairs=   8, rRMSD=0.034 ( -3.11)), URMSD= 0.049 (rURMSD=0.049)
                
                # Remove spaces after = 
                line = re.sub("= +", "=", line )
                
                # Remove trailing . after numbers
                line = re.sub("\. +", " ", line )
                
                # Now remove commas and brackets
                line = line.replace(",","")
                line = line.replace("("," ")
                
                fields = line.split()
                
                label, value = fields[0].split( "=" )
                assert label == "RMSD"
                d.rmsd = float( value )
                
                label, value = fields[1].split( "=" )
                assert label == "Pairs"
                d.pairs = int( value )
                
        return d

    def tm(self, modelName ):
        """"""
        for d in self.data:
            if d.modelName == modelName:
                return d.tm
        assert False
        return
            
    def rmsd(self, modelName ):
        """"""
        for d in self.data:
            if d.modelName == modelName:
                return d.rmsd
        assert False
        return

    def maxsubSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.maxsub, reverse=reverse )
     
    def runCompareDirectory(self, nativePdb=None, modelsDirectory=None, logfile=None ):
        
        # Generate the list of models
        pdblist = os.path.join( self.workdir, "models.list")
        with open( pdblist, 'w' ) as f:
            f.write( os.linesep.join( glob.glob( os.path.join( modelsDirectory, 'S_*.pdb' ) ) ) )
            
        # Run Maxcluster
        cmd = [ self.maxclusterExe, "-e", nativePdb, "-l", pdblist, ]
        retcode = ample_util.run_command( cmd, logfile=logfile, dolog=True )
        
        if retcode != 0:
            msg = "non-zero return code for maxcluster in runMaxcluster!"
            #logging.critical( msg )
            raise RuntimeError, msg
    
        return
     
    def tmSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.tm, reverse=reverse )

