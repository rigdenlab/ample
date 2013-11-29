#!/usr/bin/env python
'''
Created on 19 Jul 2013

@author: jmht

Data we will be collecting:

Target:
length in AA
resolution in A
secondary structure (% and per AA)
?radius of gyration

Rosetta Models:
Score
maxsub cf target

Cluster:

Ensemble:
number of models (and their scores)
truncation level
number residues
side_chain_treatment

Solution
number (and identity? of residues)

shelxe rebuild:
* CC
* av. fragment length
* RMSD to native
* Maxsub to native
* TM to native

remac refined result
* reforigin score to native
* rmsd to native
* maxsub to native
* TM to native

TO THINK ABOUT
* multiple models/chains in native
* multiple chains in solution (e.g. 3PCV)

Strategy:
for each PDB structure
Extract Data (e.g. resolution,length) from native PDB.
Extract DSSP data from native DSSP file.
Extract secondary structure prediction (SSP) from fragment prediction.
Compare DSSP and SSP data.
Extract Maxsub and RMSD data from ROSETTA score file for all models.
For each result:
Determine the ensemble that generated the result.
Extract data on ensemble from serialised file for run.
Extract data on MRBUMP MR job from logfiles

Calculate RMSD:
Standardise native PDB (most probable conformation, HETATM etc.).
Loop over each chain in the native extracting chain to file:
Loop over each chain in refined PDB extracting to file
Calculate map between native and model - NEED FULL FILE?
strip native PDB to C-alpha atoms
create a PDB containing only those ATOMS from native that are in refined PDB
and renaming the chain and renumbering the atoms to match the refined PDB.
Use reforigin to compare.





'''

import cPickle
import csv
import glob
import os
import re
import shutil
import sys
import types
import unittest

#sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
sys.path.append("/opt/ample-dev1/python")

import ample_util
import contacts
import dssp
import mrbump_results
import pdb_edit
import phaser_parser
import reforigin
import residue_map
import shelxe_log

class AmpleResult(object):
    """Results for an ample solution"""
    
    def __init__(self):



        # The attributes we will be holding
        self.orderedAttrs = [ 
                              'pdbCode',
                              'title',
                              'fastaLength',
                              'numChains',
                              'estChainsASU',
                              'resolution',
                              'solventContent',
                              'matthewsCoefficient',
                              'ss_pred',
                              'ss_pred_str',
                              'ss_dssp',
                              'ss_dssp_str',
                              'resultDir',
                              'spickerClusterSize',
                              'spickerClusterCentroid',
                              'ensembleName',
                              'ensembleNumModels',
                              'ensembleNumResidues',
                              'ensemblePercentModel',
                              'ensembleSideChainTreatment',
                              'ensembleRadiusThreshold',
                              'ensembleTruncationThreshold',
                              'ensembleNativeRmsd',
                              'ensembleNativeTM',
                              'mrProgram',
                              'phaserLLG',
                              'phaserTFZ',
                              'phaserTime',
                              'phaserKilled',
                              'molrepScore',
                              'molrepTime',
                              'reforiginRmsd',
                              'floatingOrigin',
                              'csymmatchOrigin',
                              'contactData',
                              'contactOrigin',
                              'numContacts',
                              'inregisterContacts',
                              'ooregisterContacts',
                              'backwardsContacts',
                              'goodContacts',
                              'nocatContacts',
                              'helixSequence',
                              'lenHelix',
                              'rfact',
                              'rfree',
                              'solution',
                              'shelxeCC',
                              'shelxeAvgChainLength'
                              ]
        
        # The matching titles
        self.orderedTitles = [  
                                "PDB Code",
                                "Title",
                                "Fasta Length",
                                "Number of Chains",
                                "Est. Chains in ASU",
                                "Resolution",
                                "Solvent Content",
                                "Matthews Coefficient",
                                "SS_Prediction Data",
                                "SS Prediction",
                                "DSSP Data",
                                "DSSP Assignment",
                                "Result Dir",
                                'Spicker cluster size',
                                'Spicker Centroid',
                                "Ensemble name",
                                "Ensemble num models",
                                "Ensemble num residues",
                                "Ensemble % of Model",
                                "Ensemble side chain",
                                "Ensemble radius thresh",
                                "Ensemble truncation thresh",
                                'Ensemble Native Rmsd',
                                'Ensemble Native TM',
                                "MR program",
                                "Phaser LLG",
                                "Phaser TFZ",
                                "Phaser Time",
                                "Phaser Killed",
                                "Molrep Score",
                                "Molrep Time",
                                "Reforigin RMSD",
                                "Floating Origin",
                                "Csymmatch Origin",
                                "Contact Data",
                                "Contact origin",
                                "Number of contacts",
                                "In register contacts",
                                "Out of register contacts",
                                "Backwards contacts",
                                "Good contacts",
                                "Uncategorised contacts",
                                "Helix sequence",
                                "Helix length",
                                "Rfact",
                                "Rfree",
                                "Solution",
                                "Shelxe CC",
                                "Shelxe avg. chain length"
                                 ]

        # Things not to output
        self.skip = [ "resultDir", "ss_pred", "ss_dssp", "contactData" ]
        
        # Set initial values
        for a in self.orderedAttrs:
            setattr( self, a, None )
        
        return
    
    def valuesAsList(self):
        return [ getattr(self, a) for a in self.orderedAttrs if a not in self.skip ]
    
    def titlesAsList(self):
#         l = []
#         for i, t in enumerate( self.orderedTitles ):
#             if self.orderedAttrs[i] not in self.skip:
#                 l.append( t )
#         return l
        return [ t for i, t in enumerate( self.orderedTitles ) if self.orderedAttrs[i] not in self.skip ]
    
    def asDict(self):
        """Return ourselves as a dict"""
        d = {}
        for a in self.orderedAttrs:
            d[ a ] = getattr(self, a)
        return d
    
    def __str__(self):
        
        s = ""
        for i, t in enumerate( self.orderedTitles ):
            if self.orderedAttrs[i] in self.skip:
                continue
            s += "{0:<26} : {1}\n".format( t, getattr( self, self.orderedAttrs[i] ) )
        return s
    
# End AmpleResult

class CompareModels(object):
    """Class to compare two models - currently with maxcluster"""
    
    def __init__(self, refModel, targetModel, workdir=None ):
        
        self.workdir = workdir
        
        self.refModel = refModel
        self.targetModel = targetModel
        
        self.grmsd = None
        self.maxsub = None
        self.pairs = None
        self.rmsd = None
        self.tm = None
        
        pdbedit = pdb_edit.PDBEdit()
        
        #print "CompareModels refModel: {0}  targetModel: {1}".format( refModel, targetModel )
        
        # If the rebuilt models is in multiple chains, we need to create a single chain
        nativeInfo = pdbedit.get_info( self.targetModel )
        if len( nativeInfo.models[0].chains ) > 1:
            #print "Coallescing targetModel into a single chain"
            n = os.path.splitext( os.path.basename( self.targetModel ) )[0]
            targetModel1chain = os.path.join( workdir, n+"_1chain.pdb" )
            pdbedit.to_single_chain( self.targetModel, targetModel1chain )
            self.targetModel = targetModel1chain
            
        # If the reference model is in multiple chains, we need to create a single chain
        nativeInfo = pdbedit.get_info( self.refModel )
        if len( nativeInfo.models[0].chains ) > 1:
            #print "Coallescing refModel into a single chain"
            n = os.path.splitext( os.path.basename( self.refModel ) )[0]
            refModel1chain = os.path.join( workdir, n+"_1chain.pdb" )
            pdbedit.to_single_chain( self.refModel, refModel1chain )
            self.refModel = refModel1chain
        
        self.run()
    
        return
    
    def run(self):
        
        n = os.path.splitext( os.path.basename( self.targetModel ) )[0]
        logfile = os.path.join( self.workdir, n+"_maxcluster.log" )
        
        # Run maxcluster in sequence independant mode
        cmd="/opt/maxcluster/maxcluster -in -e {0} -p {1}".format( self.targetModel, self.refModel ).split()
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False)
        
        if retcode != 0:
            raise RuntimeError,"Error running maxcluster!"
        
        self.parse_maxcluster_log( logfile )
        
        alignrsm = os.path.join( self.workdir, "align.rsm")
        
        n = os.path.splitext( os.path.basename( self.targetModel ) )[0]
        rootname = os.path.join( self.workdir, n )
        self.split_alignrsm( alignrsm=alignrsm, rootname=rootname )
        
        return
    
    def parse_maxcluster_log( self, logfile ):
        """Extract nativeInfo - assumes it completers in 1 iteration"""
        
        
        def _get_float( istr ):
            # Needed as we sometimes get spurious characters after the last digit
            nums = [ str(i) for i in range(10 ) ]
            if istr[-1] not in nums:
                istr = istr[:-1]
            return float(istr)
        
        for line in open( logfile , 'r' ):
            if line.startswith("Iter"):
                # colon after int
                iternum = int( line.split()[1][:-1] )
                if iternum > 1:
                    raise RuntimeError,"More than one iteration - no idea what that means..."
                
                if line.find( " Pairs=") != -1:
                    # e.g.: Iter 1: Pairs= 144, RMSD= 0.259, MAXSUB=0.997. Len= 144. gRMSD= 0.262, TM=0.997
                    # Need to remove spaces after = sign as the output is flakey - do it twice for safety
                    tmp = line.replace("= ","=")
                    tmp = tmp.replace("= ","=")
                    for f in tmp.split():
#                         if f.startswith("Pairs"):
#                             self.pairs = int( f.split("=")[1] )
                        if f.startswith("RMSD"):
                            self.rmsd = _get_float( f.split("=")[1] )
                        if f.startswith("MAXSUB"):
                            self.maxsub = _get_float( f.split("=")[1] )
                        if f.startswith("gRMSD"):
                            self.grmsd = _get_float( f.split("=")[1] )
                        if f.startswith("TM"):
                            self.tm = _get_float( f.split("=")[1] )
                            
                    # Bail out as we should be done
                    break
            
        return
    
    def split_alignrsm(self, alignrsm=None, rootname=None):
         
        # Order is experiment then assignment
 
        
        efile = rootname+"_experiment.pdb"
        pfile = rootname+"_prediction.pdb"
                 
        f = open( alignrsm, 'r' )
        line = f.readline()
         
        reading = False
        gotExp = False
        lines = []
        for line in open( alignrsm, 'r' ):
            
            # End of one of the files
            if line.startswith("TER"):
                lines.append( line )
                if not gotExp:
                    gotExp=True
                    f = open( efile, 'w' )
                else:
                    f = open( pfile, 'w' )
                    
                f.writelines( lines )
                f.close()
                lines = []
                reading=False
            
            if line.startswith("REMARK"):
                if not reading:
                    reading = True
                    
            if reading:
                lines.append( line )
        
        return
# End CompareModels



class EnsemblePdbParser(object):
    """
    Class to mine information from an ensemble pdb
    """

    def __init__(self,pdbfile):

        self.pdbfile = pdbfile
        
        self.centroidModelName = None
        self.modelNames = []
        self.models = []

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)

        capture=False
        for line in open(self.pdbfile, 'r'):
            
            if line.startswith( "REMARK   MODEL" ) and not capture:
                capture=True
                
            if capture and not line.startswith( "REMARK   MODEL" ):
                break
            
            if capture:
                fields = line.split()
                self.models.append( fields[3] )
                
        
        if not len( self.models ):
            raise RuntimeError,"Failed to get any models from ensemble!"
        
        for m in self.models:
            self.modelNames.append( os.path.splitext( os.path.basename( m ) )[0] )
            
        self.centroidModelName = self.modelNames[0]

        return

class MaxclusterData(object):
    
    def __init__(self):
        self.pairs = None
        self.rmsd = None
        self.maxsub = None
        self.tm = None
        self.msi = None
        self.modelName = None
        self.pdb = None
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

class MaxclusterComparator(object):
    """
        
    # Extract the first chain from the nativePdb
    
    # Create a residueSequenceMap and see if the residues match
    
    # If not use keep_matching to create a nativePdb that has the correct residue sequence`
    
    # Run Maxcluster to compare the models to the native
    """
    
    def __init__(self, nativePdb, modelsDirectory, workdir=None ):
        
        
        self.data = []
        self.modelsDirectory = modelsDirectory
        
        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
        
        self.maxclusterLogfile = os.path.join( self.workdir, "maxcluster.log" )
        self.nativePdb = self.prepareNative( nativePdb )
        
        #self.maxclusterExe = "/Users/jmht/Documents/AMPLE/programs/maxcluster"
        self.maxclusterExe = "/opt/maxcluster/maxcluster"
        
        self.runMaxcluster()
        
        self.parseLog( )
        
    def prepareNative(self, nativePdb ):
        """do stuff"""
        
        # Find out how many chains and extract the first if > 1
        
        PE = pdb_edit.PDBEdit()
        info = PE.get_info( nativePdb )
        if len( info.models ) > 1:
            raise RuntimeError,"More than 1 model."
        
        # Check if > 1 chain
        chainID=None
        if len( info.models[0].chains ) > 1:
            chainID=info.models[0].chains[0]
        
        # Standardise the native and extract the chain if > 1
        n = os.path.splitext( os.path.basename( nativePdb ) )[0]
        nativePdbStd = "{0}_std.pdb".format( n )
        PE.standardise(nativePdb, nativePdbStd, chain=chainID )
        nativePdb = nativePdbStd
        
        # 1 chain so check if the residues and sequences match
        model = os.path.join( self.modelsDirectory, "S_00000001.pdb" )
        resMap = residue_map.residueSequenceMap( nativePdb, model )
        
        #if resMap.nativeResSeq != resMap.modelResSeq or resMap.nativeSequence != resMap.model.Sequence:
        if not resMap.resSeqMatch():
            
            # We need to create a copy of the native with numbering matching the model
            n = os.path.splitext( os.path.basename( nativePdb ) )[0]
            nativeRenumber = "{0}_renumber.pdb".format( n )
            PE.match_resseq( targetPdb=nativePdb, outPdb=nativeRenumber, resMap=resMap )
            
            nativePdb = nativeRenumber
        
        return nativePdb
        
    def parseLog(self):
        
        self.data = []
        
        #INFO  : 1000. 2XOV_clean_ren.pdb vs. /media/data/shared/TM/2XOV/models/S_00000444.pdb  Pairs=  36, RMSD= 3.065, MaxSub=0.148, TM=0.192, MSI=0.148
        for line in open( self.maxclusterLogfile, 'r' ):
            
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


    def tm(self, modelName ):
        """"""
        for d in self.data:
            if d.modelName == modelName:
                return d.tm
            
    def rmsd(self, modelName ):
        """"""
        for d in self.data:
            if d.modelName == modelName:
                return d.rmsd

    def maxsubSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.maxsub, reverse=reverse )
     
    def runMaxcluster(self):
        
        # Generate the list of models
        pdblist = os.path.join( self.workdir, "models.list")
        with open( pdblist, 'w' ) as f:
            f.write( os.linesep.join( glob.glob( os.path.join( self.modelsDirectory, 'S_*.pdb' ) ) ) )
            
        # Run Maxcluster
        cmd = [ self.maxclusterExe, "-e", self.nativePdb, "-l", pdblist, ]
        retcode = ample_util.run_command( cmd, logfile=self.maxclusterLogfile, dolog=False )
        
        if retcode != 0:
            msg = "non-zero return code for maxcluster in runMaxcluster!"
            #logging.critical( msg )
            raise RuntimeError, msg
     
    def tmSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.tm, reverse=reverse )


class MolrepLogParser(object):
    """
    Class to mine information from a 
    """

    def __init__(self,logfile):

        self.logfile = logfile

        self.score=None
        self.tfScore=None
        self.wrfac=None
        self.time = None

        self.parse()

        return

    def parse(self):
        """This just drops through reading each summary and so we are left with the last one"""

        fh=open(self.logfile)

        line=fh.readline()
        while line:
            
            if "--- Summary ---" in line:
                # really scrappy - just skip 3 and take whatever comes next
                fh.readline()
                fh.readline()
                fh.readline()
                line=fh.readline()
                fields = line.split()
                if len(fields) != 14:
                    raise RuntimeError,"Error reading summary for line: {0}".format( line )
                
                self.tfScore = float( fields[10] )
                self.wrfac = float( fields[11] )
                self.score= float( fields[12] )
                
            if line.startswith( "Times: User:" ):
                fields = line.split()
                time = fields[6]
                m,s = time.split(":")
                self.time = int(m)*60 + int(s)

            line=fh.readline()
        fh.close()

        return

class MrbumpLogParser(object):
    """
    Class to mine information from a 
    """

    def __init__(self,logfile):

        self.logfile = logfile

        self.noResTarget=0
        self.noChainsTarget=0
        self.resolution=0.0

        self.parse()

        return

    def parse(self):
        """parse"""

        fh=open(self.logfile)

        line=fh.readline()
        while line:
            if "Number of residues:" in line:
                self.noResTarget=int( line.split()[-1] )
            if "Estimated number of molecules to search for in a.s.u.:" in line:
                self.noChainsTarget=int( line.split()[-1] )
            if "Resolution of collected data (angstroms):" in line:
                self.resolution=float( line.split()[-1] )
            line=fh.readline()
        fh.close()

        return

class PsipredParser(object):
    """
    Class to mine information from a psipred format file
    """

    def __init__(self,pfile):

        self.psipredfile = pfile
        
        self.residues = []
        self.assignment = []
        self.percentH = None
        self.percentC = None
        self.percentE = None

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)

        self.residues = []
        self.assignment = []
        for i, line in enumerate( open(self.psipredfile, 'r') ):
            # Assume first 2 lines are not important
            if i < 2:
                continue
            
            if not line:
                break
            
            fields = line.split()
            if len(fields) != 6:
                raise RuntimeError,"Wrong number of fields in line: ",line
            
            idx, resid, pred = int(fields[0]), fields[1], fields[2]
            if 1 == 3 and idx != 1:
                raise RuntimeError,"Got wrong index for first data: ",line
            
            self.residues.append( resid )
            self.assignment.append( pred )
            
        
        if not len(self.residues) or not len( self.assignment):
            raise RuntimeError,"Got no assignment!"
        
        nH = 0
        nC = 0
        nE = 0
        for p in self.assignment:
            if p == "H":
                nH += 1
            elif p == "C":
                nC += 1
            elif p == "E":
                nE += 1
            else:
                raise RuntimeError,"Unrecognised assignment: {0}".format( p )
            
        self.percentC = float(nC) / len(self.assignment) * 100
        self.percentH = float(nH) / len(self.assignment) * 100
        self.percentE = float(nE) / len(self.assignment) * 100
            
            
        return

    def asDict(self):
        d = {}
        d['assignment'] = self.assignment
        d['residues'] = self.residues
        d['percentC'] = self.percentC
        d['percentE'] = self.percentE
        d['percentH'] = self.percentH
        
        return d
    
class RefmacLogParser(object):
    """
    Class to mine information from a phaser log
    """

    def __init__(self,logfile):

        self.logfile = logfile

        self.initRfree=1.0
        self.finalRfree=1.0
        self.initRfactor=1.0
        self.finalRfactor=1.0
        self.noResModel=0
        self.noChainsModel=0

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)
        fh = open(self.logfile, 'r')

        CAPTURE=False

        line=fh.readline()
        while line:
            if "Number of residues :" in line:
                self.noResModel=int( line.split()[-1] )
            if "Number of chains   :" in line:
                self.noChainsModel=int( line.split()[-1] )
            if CAPTURE:
                if "R free" in line:
                    self.initRfree=float( line.split()[-2] )
                    self.finalRfree=float( line.split()[-1] )
                    CAPTURE=False
                if "R factor" in line:
                    self.initRfactor=float( line.split()[-2] )
                    self.finalRfactor=float( line.split()[-1] )
            if " $TEXT:Result: $$ Final results $$" in line:
                CAPTURE=True
            line=fh.readline()
        fh.close()

        return

class RosettaScoreData(object):
    
    def __init__(self):
        self.score = None
        self.rms = None
        self.maxsub = None
        self.description = None
        self.model = None
        return

class RosettaScoreParser(object):
    
    def __init__(self, directory ):
        
        self.directory = directory
        
        self.avgScore = None
        self.topScore = None
        self.avgRms = None
        self.topRms = None
        self.avgMaxsub = None
        self.topMaxsub = None
        
        self.data = []
        
        score_file = os.path.join( directory, "score.fsc")
        self.parseFile( score_file )
        
    def parseFile(self, score_file ):
        for i, line in enumerate( open(score_file, 'r') ):
            
            if i == 0:
                continue
    
            line = line.strip()
            if not line: # ignore blank lines - not sure why they are there...
                continue
            
            d = RosettaScoreData()
            
            fields = line.split()
            d.score = float(fields[1])
            d.rms = float(fields[26])
            d.maxsub = float(fields[27])
            d.description = fields[ 31 ]
            #pdb = fields[31]
            
            d.model = os.path.join( self.directory, d.description+".pdb" )
            
            self.data.append( d )
        
        avg = 0
        self.topScore = self.data[0].score
        for d in self.data:
            avg += d.score
            if d.score < self.topScore:
                self.topScore = d.score
        self.avgScore  = avg / len(self.data)
        
        avg = 0
        self.topRms = self.data[0].rms
        for d in self.data:
            avg += d.rms
            if d.rms < self.topRms:
                self.topRms = d.rms
        self.avgRms  = avg / len(self.data)
        
        avg = 0
        self.topMaxsub = self.data[0].maxsub
        for d in self.data:
            avg += d.maxsub
            if d.maxsub > self.topMaxsub:
                self.topMaxsub = d.maxsub
        self.avgMaxsub  = avg / len(self.data)
        
    def maxsubSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.maxsub, reverse=reverse )
     
    def rmsSorted(self, reverse=True ):
        return sorted( self.data, key=lambda data: data.rms, reverse=reverse )
    
    def rms(self, name):
        for d in self.data:
            if d.description == name:
                return d.rms
            
    def maxsub(self, name):
        for d in self.data:
            if d.description == name:
                return d.maxsub
    
    def __str__(self):
        s = "Results for: {0}\n".format(self.name)
        s += "Top score : {0}\n".format( self.topScore )
        s += "Avg score : {0}\n".format( self.avgScore )
        s += "Top rms   : {0}\n".format( self.topRms )
        s += "Avg rms   : {0}\n".format( self.avgRms )
        s += "Top maxsub: {0}\n".format( self.topMaxsub )
        s += "Avg maxsub: {0}\n".format( self.avgMaxsub )
        return s

class Test(unittest.TestCase):


    def testShelxeLogParser(self):
        logfile = "/media/data/shared/TM/2BHW/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_9.355791_rad_3_molrep_mrbump/" + \
        "data/loc0_ALL_poly_ala_trunc_9.355791_rad_3/unmod/mr/molrep/build/shelxe/shelxe_run.log"
        
        p = shelxe_log.ShelxeLogParser( logfile )
        self.assertEqual(37.26, p.CC)
        self.assertEqual(7, p.avgChainLength)
        self.assertEqual(9, p.maxChainLength)
        self.assertEqual(1, p.cycle)
        self.assertEqual(14, p.numChains)
    
    

#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()

if __name__ == "__main__":
    
    pickledResults=False
    CLUSTERNUM=0
    rundir = "/home/jmht/Documents/test/TM"
    rundir = "/home/jmht/Documents/test/CC/run2"
    #dataRoot = "/media/data/shared/TM"
    dataRoot = "/media/data/shared/coiled-coils"
    
    os.chdir( rundir )
    
    if pickledResults:
        pfile = os.path.join( dataRoot,"results.pkl" )
        with open( pfile ) as f:
            resultsDict = cPickle.load( f  )
    
    allResults = []
    
    for pdbcode in [ l.strip() for l in open( os.path.join( dataRoot, "dirs.list") ) if not l.startswith("#") ]:
    #for pdbcode in sorted( resultsDict.keys() ):
    #for pdbcode in [ "1M3W" ]:
        
        workdir = os.path.join( rundir, pdbcode )
        if not os.path.isdir( workdir ):
            os.mkdir( workdir )
        os.chdir( workdir )
            
        print "\nResults for ",pdbcode
        
        # Directory where all the data for this run live
        dataDir = os.path.join( dataRoot, pdbcode )
        
        # Get the path to the original pickle file
        pfile = os.path.join( dataDir, "ROSETTA_MR_0/resultsd.pkl")
        with open( pfile ) as f:
            ampleDict = cPickle.load( f  )
    
        # First process all stuff that's the same for each structure
        
        # Get path to native Extract all the nativeInfo from it
        nativePdb = os.path.join( dataDir, "{0}.pdb".format( pdbcode ) )
        pdbedit = pdb_edit.PDBEdit()
        nativeInfo = pdbedit.get_info( nativePdb )
        
        # First check if the native has > 1 model and extract the first if so
        if len( nativeInfo.models ) > 1:
            print "nativePdb has > 1 model - using first"
            nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=workdir )
            pdbedit.extract_model( nativePdb, nativePdb1, modelID=nativeInfo.models[0].serial )
            nativePdb = nativePdb1
            
        # Standardise the PDB to rename any non-standard AA, remove solvent etc
        nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=workdir )
        pdbedit.standardise( nativePdb, nativePdbStd )
        nativePdb = nativePdbStd
        
        # Get the new Info about the native
        nativeInfo = pdbedit.get_info( nativePdb )
        
        # Get the scores for the models - we use both the rosetta and maxcluster methods as maxcluster
        # requires a separate run to generate total RMSD
        scoreP = RosettaScoreParser( os.path.join( dataDir, "models") )
        maxComp = MaxclusterComparator( nativePdb, os.path.join( dataDir, "models")  )
        
        # Secondary Structure assignments
        #sam_file = os.path.join( dataDir, "fragments/t001_.rdb_ss2"  )
        psipred_file = os.path.join( dataDir, "fragments/t001_.psipred_ss2"  )
        psipredP = PsipredParser( psipred_file )
        dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbcode  )  )
        dsspP = dssp.DsspParser( dsspLog )

        # Get hold of a full model so we can do the mapping of residues
        refModelPdb = os.path.join( dataDir, "models/S_00000001.pdb".format( pdbcode ) )
        resSeqMap = residue_map.residueSequenceMap()
        
        modelInfo = pdbedit.get_info( refModelPdb )
        
        resSeqMap.fromInfo( refInfo=modelInfo,
                                refChainID=modelInfo.models[0].chains[0],
                                targetInfo=nativeInfo,
                                targetChainID=nativeInfo.models[0].chains[0]
                                )
        
        # Loop over each result
        if pickledResults:
            results = resultsDict[ pdbcode ]
        else:
            r = mrbump_results.ResultsSummary( os.path.join( dataDir, "ROSETTA_MR_0/MRBUMP/cluster_1") )
            r.extractResults()
            results = r.results
        
        #jtest=0
        for mrbumpResult in results:
            
            #jtest += 1
            #if jtest > 1:
            #    break
            
            print "processing result ",mrbumpResult.name
            
            ar = AmpleResult()
            allResults.append( ar )
            
            ar.pdbCode = pdbcode
            ar.title = nativeInfo.title
            ar.fastaLength = ampleDict['fasta_length']
            ar.spickerClusterSize = ampleDict['spicker_results'][ CLUSTERNUM ].cluster_size
            ar.spickerClusterCentroid = os.path.splitext( os.path.basename( ampleDict['spicker_results'][ CLUSTERNUM ].cluster_centroid ) )[0]
            ar.numChains = len( nativeInfo.models[0].chains )
            ar.resolution = nativeInfo.resolution
            ar.solventContent = nativeInfo.solventContent
            ar.matthewsCoefficient = nativeInfo.matthewsCoefficient      
            ar.ss_pred = psipredP.asDict()
            ar.ss_pred_str = "C:{0:d} | E:{1:d} | H:{2:d}".format( int(psipredP.percentC),  int(psipredP.percentE), int(psipredP.percentH) )
            ar.ss_dssp = dsspP.asDict()
            ar.ss_dssp_str = "C:{0:d} | E:{1:d} | H:{2:d}".format( int(dsspP.percentC[0]),  int(dsspP.percentE[0]), int(dsspP.percentH[0]) )
            
            # yuck...
            if mrbumpResult.solution == 'unfinished':
                s = mrbumpResult.name.split("_")
                # Different for CC and TM cases
                #ensembleName = "_".join( s[2:-2] )
                ensembleName = "_".join( s[2:-1] )
            else:
                # MRBUMP Results have loc0_ALL_ prepended and  _UNMOD appended
                ensembleName = mrbumpResult.name[9:-6]
            ar.ensembleName = ensembleName
            
            #if ensembleName != "All_atom_trunc_9.06527_rad_3":
            #    continue
            
            # Extract information on the models and ensembles
            eresults = ampleDict['ensemble_results']
            got=False
            for e in ampleDict[ 'ensemble_results' ][ CLUSTERNUM ]:
                if e.name == ensembleName:
                    got=True
                    break 
        
            if not got:
                raise RuntimeError,"Failed to get ensemble results"
        
            ar.ensembleNumModels =  e.num_models
            ar.ensembleNumResidues =  e.num_residues
            ar.ensembleSideChainTreatment = e.side_chain_treatment
            ar.ensembleRadiusThreshold = e.radius_threshold
            ar.ensembleTruncationThreshold =  e.truncation_threshold
            ar.ensemblePercentModel = int( ( float( ar.ensembleNumResidues ) / float( ar.fastaLength ) ) * 100 )
            
            # Get the data on the models in the ensemble
            ensembleFile = os.path.join( dataDir, "ROSETTA_MR_0/ensembles_1", ensembleName+".pdb" )
            eP = EnsemblePdbParser( ensembleFile )
            
            ar.ensembleNativeRmsd = scoreP.rms( eP.centroidModelName )
            ar.ensembleNativeTM = maxComp.tm( eP.centroidModelName )
    
            ar.solution =  mrbumpResult.solution
            
            # No results so move on
            if mrbumpResult.solution == 'unfinished':
                ar.solution = mrbumpResult.solution
                continue
    
            # Need to remove last component as we recored the refmac directory
            resultDir = os.sep.join( mrbumpResult.resultDir.split(os.sep)[:-1] )
            ar.resultDir = resultDir
            
            ar.rfact =  mrbumpResult.rfact
            ar.rfree =  mrbumpResult.rfree
            ar.mrProgram =  mrbumpResult.program
            
            #mrbumpLog = os.path.join( dataDir, "ROSETTA_MR_0/MRBUMP/cluster_1/", "{0}_{1}.sub.log".format( ensembleName, mrbumpResult.program )  )
            mrbumpLog = os.path.join( dataDir, "ROSETTA_MR_0/MRBUMP/cluster_1/", "{0}.sub.log".format( ensembleName )  )
            mrbumpP = MrbumpLogParser( mrbumpLog )
            ar.estChainsASU = mrbumpP.noChainsTarget
            
            if mrbumpResult.program == "phaser":
                
                phaserPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
                if not os.path.isfile( phaserPdb ):
                    continue
                
                phaserP = phaser_parser.PhaserPdbParser( phaserPdb )
                ar.phaserLLG = phaserP.LLG
                ar.phaserTFZ = phaserP.TFZ
                
                phaserLog = os.path.join( resultDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(mrbumpResult.program, ensembleName) )
                phaserP = phaser_parser.PhaserLogParser( phaserLog )
                ar.phaserTime = phaserP.time
                
                placedPdb = phaserPdb
                
            else:
                molrepLog = os.path.join( resultDir, "molrep.log" )
                molrepP = MolrepLogParser( molrepLog )
                ar.molrepScore = molrepP.score
                ar.molrepTime = molrepP.time
                
                placedPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
                if not os.path.isfile( placedPdb ):
                    continue
                
            # Get the reforigin RMSD of the phaser placed model as refined with refmac
            #refinedPdb = os.path.join( resultDir, "refine", "refmac_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
            
            # debug - copy into work directory as reforigin struggles with long pathnames
            shutil.copy(placedPdb, os.path.join( workdir, os.path.basename( placedPdb ) ) )
            
            placedInfo = pdbedit.get_info( placedPdb )
            
            # Get reforigin info
            try:
                rmsder = reforigin.ReforiginRmsd( nativePdb=nativePdb,
                                                  nativePdbInfo=nativeInfo,
                                                  placedPdb=placedPdb,
                                                  placedPdbInfo=placedInfo,
                                                  refModelPdb=refModelPdb,
                                                  cAlphaOnly=True )
                ar.reforiginRmsd = rmsder.rmsd
            except Exception, e:
                print "ERROR: ReforiginRmsd with: {0} {1}".format( nativePdb, placedPdb )
                print "{0}".format( e )
                ar.reforiginRmsd = 9999
                
            #
            # SHELXE PROCESSING
            #
            # Now read the shelxe log to see how we did
            shelxeLog = os.path.join( resultDir, "build/shelxe/shelxe_run.log" )
            origShelxePdb = os.path.join( resultDir, "build/shelxe", "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
            shelxePdb = None
            if os.path.isfile( origShelxePdb):
                shelxePdb = os.path.join(workdir, os.path.basename( origShelxePdb ) )
                shutil.copy( origShelxePdb, shelxePdb   )
                
            if os.path.isfile( shelxeLog ):
                shelxeP = shelxe_log.ShelxeLogParser( shelxeLog )
                ar.shelxeCC = shelxeP.CC
                ar.shelxeAvgChainLength = shelxeP.avgChainLength
                
            # Now calculate contacts
            ccalc = contacts.Contacts()
            try:
                ccalc.getContacts( nativePdb=nativePdb,
                                   placedPdb=placedPdb,
                                   resSeqMap=resSeqMap,
                                   nativeInfo=nativeInfo,
                                   shelxePdb=shelxePdb,
                                   workdir=workdir,
                                   dsspLog=dsspLog
                                )
            except Exception, e:
                print "ERROR WITH CONTACTS: {0}".format( e )
       
            if ccalc.best:
                ar.contactData = ccalc.best
                ar.numContacts = ccalc.best.numContacts
                ar.inregisterContacts = ccalc.best.inregister
                ar.ooregisterContacts = ccalc.best.ooregister
                ar.backwardsContacts = ccalc.best.backwards
                ar.contactOrigin = ccalc.best.origin
                ar.goodContacts = ar.inregisterContacts + ar.ooregisterContacts
                ar.nocatContacts = ar.numContacts - ar.goodContacts
                ar.helixSequence = ccalc.best.helix
                if ccalc.best.helix:
                    ar.lenHelix = len( ccalc.best.helix )
                
                hfile = os.path.join( workdir, "{0}.helix".format( ensembleName ) )
                if not ccalc.writeHelixFile( hfile ):
                    print "NO HELIX FILE"
                        
                # Just for debugging
                if ar.shelxeCC >= 25 and ar.shelxeAvgChainLength >= 10:
                    # Show origin stats
                    oc = sorted(ccalc.originCompare.items(), key=lambda x: x[1], reverse=True )
                    duff=False
                    if len(oc) > 1:
                        if oc[0][1] == oc[1][1]:
                            if len(oc) > 2:
                                if oc[2][1] >= oc[1][1]*.5:
                                    duff=True
                        else:
                            if oc[1][1] >= oc[0][1]*.5:
                                duff=True
                        if duff:
                            print "OTHER ORIGINMATCHES ARE > 50%"
                            print "originCompare: ", oc
                 
                #print ar

    # End loop over results
    
    pfile = os.path.join( rundir, "ar_results.pkl")
    f = open( pfile, 'w' )
    ampleDict = cPickle.dump( allResults, f  )
    
    cpath = os.path.join( rundir, 'results.csv' )
    csvfile =  open( cpath, 'wb')
    csvwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    header=False
    for r in allResults:
        if not header:
            csvwriter.writerow( r.titlesAsList() )
            header=True
        csvwriter.writerow( r.valuesAsList() )
        
    csvfile.close()
