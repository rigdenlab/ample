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
import os
import re
import shutil
import sys
import unittest

#sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
sys.path.append("/opt/ample-dev1/python")
import ample_util
import mrbump_results
import pdb_edit



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
                              'ensembleName',
                              'ensembleNumModels',
                              'ensembleNumResidues',
                              'ensemblePercentModel',
                              'ensembleSideChainTreatment',
                              'ensembleRadiusThreshold',
                              'ensembleTruncationThreshold',
                              'ensembleNativeRmsd',
                              'ensembleNativeMaxsub',
                              'mrProgram',
                              'phaserLLG',
                              'phaserTFZ',
                              'phaserTime',
                              'molrepScore',
                              'molrepTime',
                              'reforiginRmsd',
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
                                "Ensemble name",
                                "Ensemble num models",
                                "Ensemble num residues",
                                "Ensemble % of Model",
                                "Ensemble side chain",
                                "Ensemble radius thresh",
                                "Ensemble truncation thresh",
                                'Ensemble Native Rmsd',
                                'Ensemble Native Maxsub',
                                "MR program",
                                "Phaser LLG",
                                "Phaser TFZ",
                                "Phaser Time",
                                "Molrep Score",
                                "Molrep Time",
                                "Reforigin RMSD",
                                "Rfact",
                                "Rfree",
                                "Solution",
                                "Shelxe CC",
                                "Shelxe avg. chain length"
                                 ]

        # Things not to output
        self.skip = [ "resultDir", "ss_pred", "ss_dssp" ]
        
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

class DsspParser(object):
    """
    Class 
    """

    def __init__(self,pfile):

        self.dsspfile = pfile
        
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
        
        capture=False
        for line in open(self.dsspfile, 'r'):
            
            if "#  RESIDUE" in line:
                capture=True
                continue
                
            if capture:
                #print "\"{0}\"".format(line)
                #idx = int( line[0:5].strip() )
                #resIdx = int( line[5:10].strip() )
                #chainId = line[10:12].strip()
                resName = line[12:14].strip()
                assign = line[14:17].strip()
                #print "\"{0}\"".format(line[14:17])
                
                # Only capture first chain
                if "!" in assign:
                    capture=False
                    break
                 
                self.residues.append( resName )
                self.assignment.append( assign )
                
        if not len(self.residues) or not len( self.assignment):
            raise RuntimeError,"Got no assignment!"
         
        nH = 0
        nC = 0
        nE = 0
        for p in self.assignment:
            if p == "H":
                nH += 1
            elif p == "E":
                nE += 1
            # Just assume everything else is a coil
            else:
                nC += 1
             
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

class PhaserLogParser(object):
    """
    Class to mine information from a phaser log
    """

    def __init__(self,logfile):

        self.logfile = logfile
        self.phaserLLG = None
        self.phaserTFZ = None
        self.phaserTime = None

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)
        fh = open(self.logfile, 'r')
        
        #print "Checking logfile ",self.logfile

        for line in reversed(fh.readlines()):
            if "CPU Time" in line:
                self.phaserTime=float( line.split()[-2] )
                break
        fh.close()

        fh = open(self.logfile, 'r')
        CAPTURE = False
        solline = ""
        self.phaserLLG = 0.0
        self.phaserTFZ = 0.0
        line = fh.readline()
        while line:
            if CAPTURE:
                if "SOLU SPAC" in line:
                    CAPTURE = False
                else:
                    solline += line.strip() + " "
            if  "Solution #1 annotation (history):" in line or  "Solution annotation (history):" in line:
                CAPTURE = True
            line = fh.readline()
        fh.close()

        llist = solline.split()
        llist.reverse()
        for i in llist:
            if "TFZ==" in i and "*" not in i:
                self.phaserTFZ = float(i.replace("TFZ==", ""))
                break
            if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                self.phaserTFZ = float(i.replace("TFZ=", ""))
                break

        for i in llist:
            if "LLG==" in i:
                self.phaserLLG = float(i.replace("LLG==", ""))
                break
            if "LLG=" in i and "LLG==" not in i:
                self.phaserLLG = float(i.replace("LLG=", ""))
                break
        return

class PhaserPdbParser(object):
    """
    Class to mine information from a phaser pdb file
    """

    def __init__(self,pdbfile):

        self.pdbfile = pdbfile
        self.phaserLLG = None
        self.phaserTFZ = None

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)

        for line in open(self.pdbfile, 'r'):
            if "TFZ==" in line:
                llist = line.split()
                llist.reverse()
                for i in llist:
                    if "TFZ==" in i and "*" not in i:
                        self.phaserTFZ = float(i.replace("TFZ==", ""))
                        break
                    if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                        self.phaserTFZ = float(i.replace("TFZ=", ""))
                        break
        
                for i in llist:
                    if "LLG==" in i:
                        self.phaserLLG = float(i.replace("LLG==", ""))
                        break
                    if "LLG=" in i and "LLG==" not in i:
                        self.phaserLLG = float(i.replace("LLG=", ""))
                        break
        return

class PsipredParser(object):
    """
    Class to mine information from a phaser pdb file
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

class ReforiginRmsd(object):
    """Class to use reforigin to determine how well the model was placed.
    """
    
    def __init__( self, nativePdb, refinedPdb, refModelPdb):
        
        
        self.cAlphaOnly = True # Whether to only compare c-alpha atoms
        self.rmsd = None
        self.bestNativeChain = None
        self.bestRefinedChain = None
        self.bestReforiginPdb = None
        self.refModelPdb = refModelPdb
        
        self.run( nativePdb, refinedPdb )
        
    def calc_reforigin_rmsd( self, refpdb=None, targetpdb=None, outpdb=None, DMAX=100 ):
        
        assert refpdb and targetpdb and outpdb
        
        # HACK - REFORIGIN has a limit on the length of the command line, so we need to create copies of inputfile
        # as this has the potentially longest path
        tmptarget = ample_util.tmpFileName()+".pdb"
        shutil.copy(targetpdb, tmptarget)
        
        logfile = outpdb+".log"
        cmd="reforigin xyzin {0} xyzref {1} xyzout {2} DMAX {3}".format( tmptarget, refpdb, outpdb, DMAX ).split()
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        else:
            os.unlink( tmptarget )
        
        # Parse logfile to get RMSD
        rms = None
        for line in open( logfile, 'r' ):
            if line.startswith("RMS deviation:"):
                rms = float( line.split()[-1] )
                break
        
        if not rms:
            raise RuntimeError, "Error extracting RMS from logfile: {0}".format( logfile )
        
        return rms

    def reforigin_rmsd( self, refinedPdb, nativePdb, nativeChainID=None ):
        """Use reforigin to calculate rmsd between native and refined"""
        
        workdir=os.getcwd()
        
        # Calculate the RefSeqMap - need to do this before we reduce to c-alphas
        PE = pdb_edit.PDBEdit()
        resSeqMap = PE.get_resseq_map( nativePdb, self.refModelPdb )
        
        # Find out if there are extra atoms in the model that we need to remove
        extra = resSeqMap.modelExtra()
        if len(extra):
            
            #print "GOT EXTRA ",extra
            
            n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
            refinedPdbCut = os.path.join( workdir, n+"_cut.pdb" )
            
            logfile = "{0}.log".format( refinedPdb )
            cmd="pdbcur xyzin {0} xyzout {1}".format( refinedPdb, refinedPdbCut ).split()
            
            # Build up stdin - I'm too thick to work out the selection syntax for a discrete list
            stdin = ""
            for e in extra:
                stdin += "delresidue {0}\n".format( e )
            
            retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=workdir, dolog=False, stdin=stdin)
            
            if retcode == 0:
                # remove temporary files
                os.unlink(logfile)
            else:
                raise RuntimeError,"Error deleting residues {0}".format( extra )
            
            refinedPdb = refinedPdbCut
            
        
        if self.cAlphaOnly:
            # If only alpha atoms are required, we create a copy of the model with only alpha atoms
            n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
            tmp = os.path.join( workdir, n+"_cAlphaOnly.pdb" )
            PE.calpha_only( refinedPdb, tmp )
            refinedPdb = tmp
        else:
            # Strip down to backbone atoms
            n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
            tmp = os.path.join( workdir, n+"_backbone.pdb" )
            PE.backbone( refinedPdb, tmp  )
            refinedPdb = tmp

        # Now create a PDB with the matching atoms from native that are in refined
        n = os.path.splitext( os.path.basename( nativePdb ) )[0]
        nativePdbMatch = os.path.join( workdir, n+"_matched.pdb" )
        PE.keep_matching( refpdb=refinedPdb, targetpdb=nativePdb, outpdb=nativePdbMatch, resSeqMap=resSeqMap )
        
        # Now get the rmsd
        n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
        reforiginOut = os.path.join( workdir, n+"_chain{0}_reforigin.pdb".format( nativeChainID ) )
        rms = self.calc_reforigin_rmsd( refpdb=nativePdbMatch, targetpdb=refinedPdb, outpdb=reforiginOut )
        return ( rms, reforiginOut )
    
    def run( self, nativePdb, refinedPdb ):
        """For now just save lowest rmsd - can look at collecting more nativeInfo later
        
        Currently we assume we are only given one model and that it has already been standardised.
        """

        # Run a pass to find the # chains
        pdbedit = pdb_edit.PDBEdit()
        refinedInfo = pdbedit.get_info( refinedPdb )
        nativeInfo = pdbedit.get_info( nativePdb )
        native_chains = nativeInfo.models[ 0 ].chains
        refined_chains = refinedInfo.models[ 0 ].chains # only ever one model in the refined pdb
        
        #print "got native chains ", native_chains
        #print "got refined chains ", refined_chains
            
        rmsds = {} # dict of rmsd -> ( chainIDnative, chainIDrefined, reforiginLogfile )
        
        # Match each chain in native against refined and pick the best
        for nativeChainID in native_chains:
            
            #print "native_chain: {0}".format( nativeChainID )
                    
            if len( native_chains ) == 1:
                # Don't need to do owt as we are just using the native as is
                nativeChainPdb = nativePdb
            else:
                
                # Extract the chain from the pdb
                n = os.path.splitext( os.path.basename( nativePdb ) )[0]
                nativeChainPdb = os.path.join( workdir, n+"_chain{0}.pdb".format( nativeChainID ) ) 
                pdbedit.extract_chain( nativePdb, nativeChainPdb, chainID=nativeChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
            
            for refinedChainID in refined_chains:
                
                #print "refined_chain: {0}".format( refinedChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
                
                # Extract the chain from the pdb
                n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
                refinedChainPdb = os.path.join( workdir, n+"_chain{0}.pdb".format( refinedChainID ) ) 
                pdbedit.extract_chain( refinedPdb, refinedChainPdb, chainID=refinedChainID, newChainID=nativeChainID, cAlphaOnly=self.cAlphaOnly )
                
                #print "calculating for {0} vs. {1}".format( refinedChainID, nativeChainID  )
                #print "calculating for {0} vs. {1}".format( refinedChainPdb, nativeChainPdb  )
                rmsd, refPdb  = self.reforigin_rmsd( refinedChainPdb, nativeChainPdb, nativeChainID=nativeChainID )
#                 try:
#                     rmsd, logfile  = self.reforigin_rmsd( refinedChainPdb, nativeChainPdb, nativeChainID=nativeChainID )
#                 except:
#                     print "GOT REFORIGIN ERROR for {0},{1},{2}".format( refinedChainPdb, nativeChainPdb, nativeChainID )
#                     rmsd = 99999
#                     logfile= None
                #print "got rmsd chain ",rmsd
                
                rmsds[ rmsd ] = ( nativeChainID, refinedChainID, refPdb )
                
        # End loop over chains
        # Now pick the best...
        rmsd = sorted( rmsds.keys() )[ 0 ]
        #print "Got rmsds over chains: {0}".format( rmsds )
        
        self.rmsd = rmsd
        self.bestNativeChain = rmsds[ rmsd ][0]
        self.bestRefinedChain = rmsds[ rmsd ][1]
        self.bestReforiginPdb = rmsds[ rmsd ][2]
        #print "best chain rmsd is {0} for nativeChain {1} vs refinedChain {2}".format( self.rmsd, self.bestChains[0], self.bestChains[1] )
            
        return

class RosettaScoreParser(object):
    """
    Class to mine information from a rosetta score file
    """

    def __init__(self, scorefile):

        self.scorefile = scorefile
        
        self.d = {}

        self.parse()

        return

    def parse(self):
        """parse"""

        for i, line in enumerate( open( self.scorefile ) ):
            if i==0:
                continue 
            fields = line.split()
            rmsd = float( fields[ 26 ] )
            maxsub = float( fields[ 27 ] )
            description = fields[ 31 ]
            self.d[ description ] = ( rmsd, maxsub )

        return

class ShelxeLogParser(object):
    """
    Class to mine information from a shelxe log
    """
    
    def __init__(self,logfile):
        
        self.logfile = logfile
        self.CC = None
        self.avgChainLength = None
        
        self.parse()
        
        return
        
    def parse(self):
        """Parse a shelxe log file to get the CC and average Chain length
        """
        
        cycleData = [] # List (CC,avgChainLength) tuples - ordered by cycle
        fh = open( self.logfile, 'r')
        
        line = fh.readline()
        while line:
            
            # find should be quicker then re match
            if line.find("residues left after pruning, divided into chains as follows:") != -1:
                (cc, avgChainLength) = self._parseCycle(fh)
                cycleData.append( (cc, avgChainLength) )
            
            
            if  line.find( "Best trace (cycle" ) != -1:
                # Expecting:
                #  "Best trace (cycle   1 with CC 37.26%) was saved as shelxe-input.pdb"
                cycle = int( re.search("\s\d+\s", line).group(0) )
                cc = float( re.search("\s\d+\.\d+", line).group(0) )
                
                # Check it matches
                if cycleData[ cycle-1 ][0] != cc:
                    raise RuntimeError,"Error getting final CC!"
                
                self.CC =  cycleData[ cycle-1 ][0]
                self.avgChainLength = cycleData[ cycle-1 ][1]

            line = fh.readline()
        #End while
        
        fh.close()
        
        return
        
    def _parseCycle(self, fh):
        """
        Working on assumption each cycle contains something like the below:
<log>
           223 residues left after pruning, divided into chains as follows:
 A:   6   B:   7   C:   6   D:   6   E:   8   F:   7   G:  12   H:  12   I:   5
 J:  10   K:   6   L:   6   M:   6   N:   7   O:   6   P:   7   Q:   8   R:   6
 S:   5   T:   6   U:  10   V:   9   W:  12   X:  11   Y:   8   Z:   6   Z:   6
 Z:   6   Z:   7   Z:   6

 CC for partial structure against native data =  30.92 %
 </log>
 """
        
        lengths = []
        while True:
            
            line = fh.readline().strip()
            line = line.rstrip(os.linesep)
            if not line:
                # End of reading lengths
                break
            
            # Loop through integers & add to list
            for m in re.finditer("\s\d+", line):
                lengths.append( int(m.group(0)) )
                
        # Now calculate average chain length
        if not len( lengths ):
            raise RuntimeError, "Failed to read any fragment lengths"
        
        # Average chain lengths
        avgChainLength = sum(lengths) / int( len(lengths) )        
        
        # Here should have read the  lengths so now just get the CC
        count=0
        while True:
            line = fh.readline().strip()
            if line.startswith("CC for partial structure against native data"):
                break
            else:
                count += 1
                if count > 5:
                    raise RuntimeError,"Error parsing CC score"
            
        cc = float( re.search("\d+\.\d+", line).group(0) )
        
        return ( cc, avgChainLength )

#END ShelxeLogParser

class Test(unittest.TestCase):


    def testShelxeLogParser(self):
        logfile = "/media/data/shared/TM/2BHW/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_9.355791_rad_3_molrep_mrbump/" + \
        "data/loc0_ALL_poly_ala_trunc_9.355791_rad_3/unmod/mr/molrep/build/shelxe/shelxe_run.log"
        
        p = ShelxeLogParser( logfile )
        self.assertEqual(37.26, p.CC)
        self.assertEqual(7, p.avgChainLength)
        

#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()


if __name__ == "__main__":
    pfile = "/media/data/shared/TM/results.pkl"
    f = open( pfile )
    resultsDict = cPickle.load( f  )
    f.close()
    
    rundir = "/home/jmht/Documents/test/new"
    TMdir = "/media/data/shared/TM"
    os.chdir( rundir )
    
    allResults = []
    
    #for pdbcode in [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]:
    # fails 2UUI, 3OUF, 3PCV, 3RLB, 3U2F
    
    for pdbcode in sorted( resultsDict.keys() ):
    #for pdbcode in [ "3PCV" ]:
        
        workdir = os.path.join( rundir, pdbcode )
        if not os.path.isdir( workdir ):
            os.mkdir( workdir )
        os.chdir( workdir )
            
        print "\nResults for ",pdbcode
        
        # Directory where all the data for this run live
        datadir = os.path.join( TMdir, pdbcode )
        
        # Get the path to the original pickle file
        pfile = os.path.join( datadir, "ROSETTA_MR_0/resultsd.pkl")
        f = open( pfile )
        ampleDict = cPickle.load( f  )
        f.close()    
    
        # First process all stuff that's the same for each structure
        
        # Get path to native Extract all the nativeInfo from it
        nativePdb = os.path.join( datadir, "{0}.pdb".format( pdbcode ) )
        pdbedit = pdb_edit.PDBEdit()
        nativeInfo = pdbedit.get_info( nativePdb )
        
        # First check if the native has > 1 model and extract the first if so
        if len( nativeInfo.models ) > 1:
            print "nativePdb has > 1 model - using first"
            n = os.path.splitext( os.path.basename( nativePdb ) )[0]
            nativePdb1 = os.path.join( workdir, n+"_model1.pdb" )
            pdbedit.extract_model( nativePdb, nativePdb1, modelID=nativeInfo.models[0].serial )
            nativePdb = nativePdb1
            
        # Standardise the PDB to rename any non-standard AA, remove solvent etc
        n = os.path.splitext( os.path.basename( nativePdb ) )[0]
        nativePdbStd = os.path.join( workdir, n+"_std.pdb" )
        pdbedit.standardise( nativePdb, nativePdbStd )
        nativePdb = nativePdbStd
        
        # Secondary Structure assignments
        sam_file = os.path.join( datadir, "fragments/t001_.rdb_ss2"  )
        psipredP = PsipredParser( sam_file )
        dssp_file = os.path.join( datadir, "{0}.dssp".format( pdbcode.lower()  )  )
        dsspP = DsspParser( dssp_file )
    
        # Get the scores for the models
        scoreFile = os.path.join( datadir, "models", "score.fsc" )
        scoreP = RosettaScoreParser( scoreFile )
        
        # Loop over each result
        for mrbumpResult in resultsDict[ pdbcode ]:
            
            #print "processing result ",mrbumpResult
            
            ar = AmpleResult()
            allResults.append( ar )
            
            ar.pdbCode = pdbcode
            ar.title = nativeInfo.title
            ar.fastaLength = ampleDict['fasta_length']
            ar.numChains = len( nativeInfo.models[0].chains )
            ar.resolution = nativeInfo.resolution
            ar.solventContent = nativeInfo.solventContent
            ar.matthewsCoefficient = nativeInfo.matthewsCoefficient      
            ar.ss_pred = psipredP.asDict()
            ar.ss_pred_str = "C:{0:d} | E:{1:d} | H:{2:d}".format( int(psipredP.percentC),  int(psipredP.percentE), int(psipredP.percentH) )
            ar.ss_dssp = dsspP.asDict()
            ar.ss_dssp_str = "C:{0:d} | E:{1:d} | H:{2:d}".format( int(dsspP.percentC),  int(dsspP.percentE), int(dsspP.percentH) )
            
            # yuck...
            if mrbumpResult.solution == 'unfinished':
                s = mrbumpResult.name.split("_")
                ensembleName = "_".join( s[2:-2] )
            else:
                # MRBUMP Results have loc0_ALL_ prepended and  _UNMOD appended
                ensembleName = mrbumpResult.name[9:-6]
            ar.ensembleName = ensembleName
            
            # Extract information on the models and ensembles
            eresults = ampleDict['ensemble_results']
            got=False
            clusterNum=0
            for e in ampleDict[ 'ensemble_results' ][ clusterNum ]:
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
            ensembleFile = os.path.join( datadir, "ROSETTA_MR_0/ensembles_1", ensembleName+".pdb" )
            eP = EnsemblePdbParser( ensembleFile )
            #scoreFile = os.path.join( datadir, "models", "score.fsc" )
            #scoreP = RosettaScoreParser( scoreFile )
            ar.ensembleNativeRmsd = scoreP.d[ eP.centroidModelName ][0]
            ar.ensembleNativeMaxsub = scoreP.d[ eP.centroidModelName ][1]
    
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
            
            mrbumpLog = os.path.join( datadir, "ROSETTA_MR_0/MRBUMP/cluster_1/", "{0}_{1}.sub.log".format( ensembleName, mrbumpResult.program )  )
            mrbumpP = MrbumpLogParser( mrbumpLog )
            ar.estChainsASU = mrbumpP.noChainsTarget
            
            # Get the reforigin RMSD of the phaser placed model as refined with refmac
            refinedPdb = os.path.join( resultDir, "refine", "refmac_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
            # debug - copy into work directory as reforigin struggles with long pathnames
            if not os.path.isfile( refinedPdb ):
                # If the file is missing either phaser failed or something went horribly wrong
                continue
            shutil.copy(refinedPdb, os.path.join( workdir, os.path.basename( refinedPdb ) ) )
        
            # Get hold of a full model so we can do the mapping of residues
            refModelPdb = os.path.join( datadir, "models/S_00000001.pdb".format( pdbcode ) )
            rmsder = ReforiginRmsd( nativePdb, refinedPdb, refModelPdb )
            ar.reforiginRmsd =  rmsder.rmsd
            
            if mrbumpResult.program == "phaser":
                phaserLog = os.path.join( resultDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(mrbumpResult.program, ensembleName) )
                phaserP = PhaserLogParser( phaserLog )
                ar.phaserTime = phaserP.phaserTime
                
                phaserPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
                phaserP = PhaserPdbParser( phaserPdb )
                ar.phaserLLG = phaserP.phaserLLG
                ar.phaserTFZ = phaserP.phaserTFZ
            else:
                molrepLog = os.path.join( resultDir, "molrep.log" )
                molrepP = MolrepLogParser( molrepLog )
                ar.molrepScore = molrepP.score
                ar.molrepTime = molrepP.time
     
            # Now read the shelxe log to see how we did
            shelxeLog = os.path.join( resultDir, "build/shelxe/shelxe_run.log" )
            shelxeP = ShelxeLogParser( shelxeLog )
            ar.shelxeCC = shelxeP.CC
            ar.shelxeAvgChainLength = shelxeP.avgChainLength
            
    #         # Finally use maxcluster to compare the shelxe model with the native
    #         if False:
    #             shelxeModel = os.path.join( resultDir, "build/shelxe", "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
    #             mrbumpResult.shelxModel = shelxeModel
    #             m = CompareModels( nativePdb, shelxeModel, workdir=workdir  )
    #             mrbumpResult.shelxGrmsd = m.grmsd
    #             mrbumpResult.shelxTM = m.tm
            
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
