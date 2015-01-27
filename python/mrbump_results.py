#!/usr/bin/env ccp4-python

import copy
import cPickle
import glob
import logging
import os
import sys
import types
import unittest

# Hack to make sure we can find the modules we need
if __name__ == "__main__":
    root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
    sys.path.insert( 0, os.path.join( root, "scripts" ) )

# Our imports
import parse_arpwarp
import parse_buccaneer
import parse_phaser
#import parse_shelxe
import printTable

# We need a null logger so that we can be used without requiring a logger
class NullHandler(logging.Handler):
    def emit(self, record):
        pass
    
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

class MrBumpResult(object):
    """
    Class to hold the result of running a MRBUMP job
    """
    def __init__(self):
        """

        """
        self.jobDir = None # directory jobs ran in
        self.mrDir = None # where the actual results are
        self.pdb = None # The refmac-refined pdb
        self.name = None # The MRBUMP name
        self.ensembleName = None
        self.program = None
        self.solution = None

        self.rfact = 1.0
        self.rfree = 1.0
        self.buccRfact = 1.0 # Bucc_final_Rfact
        self.buccRfree = 1.0
        self.arpWarpRfact = 1.0 # ARP_final_Rfact/Rfree
        self.arpWarpRfree = 1.0
        self.shelxeCC = -1.0
        self.shelxeACL = -1.0
        
        # Buccaneer & arpwarp scores after final rebuild
        self.buccFinalRfact = 1.0
        self.buccFinalRfree  = 1.0
        self.arpWarpFinalRfact = 1.0
        self.arpWarpFinalRfree  = 1.0
        
        # pdb files
        self.phaserPdb=None
        self.molrepPdb=None
        self.refmacPdb=None
        self.shelxePdb=None
        self.buccaneerPdb=None
        self.arpWarpPdb=None

        # misc data
        self.phaserTFZ = None
        self.phaserLLG = None
        self.phaserTime = None
        self.phaserKilled = None
        self.molrepScore = None
        self.molrepTime = None

        self.header = [] # The header format for this table
        
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

class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """

    def __init__(self):

        self.results = []
        # List of all the possible column titles and their result object attributes
        self.title2attr = {
                            'Model_Name'       : 'name',
                            'MR_Program'       : 'program',
                            'Solution_Type'    : 'solution',
                            'final_Rfact'      : 'rfact',
                            'final_Rfree'      : 'rfree',
                            'SHELXE_CC'        : 'shelxeCC',
                            'SHELXE_ACL'       : 'shelxeACL',
                            'Bucc_final_Rfact' : 'buccRfact',
                            'Bucc_final_Rfree' : 'buccRfree',
                            'ARP_final_Rfact'  : 'arpWarpRfact',
                            'ARP_final_Rfree'  : 'arpWarpRfree',
                             }
        
        # Same as titles but Model_Name swapped for Ensemble NAme
        self.header= ['Ensemble_Name' ,
                      'MR_Program',
                      'Solution_Type',
                      'final_Rfact',
                      'final_Rfree',
                      'SHELXE_CC',
                      'SHELXE_ACL',
                      'Bucc_final_Rfact',
                      'Bucc_final_Rfree',
                      'ARP_final_Rfact',
                      'ARP_final_Rfree',
                      ]

        self.logger = logging.getLogger()
        # Add Null logger so we can be used without requiring a logger
        self.logger.addHandler(NullHandler())

        return

    def _addFailedResults(self, mrbumpDir, failed):
        """Add failures to self.results

        Args:
        failed: dict of {ensemble : result}
        header: list with header for results table
        """
        count=0
        for ensemble, reason in failed.iteritems():
            d = self.createDict()
            # name hard-coded
            d['Name'] = "loc0_ALL_" + ensemble + "_UNMOD"
            d['SearchModel_name']=ensemble 
            d['Job_dir'] = os.path.join( mrbumpDir, 'search_'+ensemble+'_mrbump' )
            d['Solution_Type'] = reason
            self.results.append( d )
            count += 1
        self.logger.debug("Added {0} MRBUMP result failures".format(count) )
        return

    def createDict(self):
        d = {}
        d['Job_dir'] = None
        d['MR_dir'] = None
        d['SearchModel_name'] = None
        d['MR_program'] = None
        d['Solution_Type'] = None
        d['final_Rfact'] = None
        d['final_Rfree'] = None
        d['SHELXE_CC'] = None
        d['SHELXE_ACL'] = None
        d['Bucc_final_Rfact'] = None
        d['Bucc_final_Rfree'] = None
        d['ARP_final_Rfact'] = None
        d['ARP_final_Rfree'] = None
        
        d['PHASER_pdb'] = None
        d['MOLREP_pdb'] = None
        d['REFMAC_pdb'] = None
        d['SHELXE_pdb'] = None
        d['Bucc_final_pdb'] = None
        d['ARP_final_pdb'] = None
        
        d['PHASER_TFZ'] = None
        d['PHASER_LLG'] = None
        d['PHASER_time'] = None
        d['PHASER_killed'] = None
        
        d['MOLREP_score'] = None
        d['MOLREP_time'] = None
        
        return d

    def extractResults( self, mrbumpDir ):
        """
        Find the results from running MRBUMP and sort them
        """
        
        mrbumpDir = os.path.abspath( mrbumpDir )
        
        if not os.path.isdir(mrbumpDir):
            self.logger.warn("extractResults - is not a valid directory: {0}".format( mrbumpDir ) )
                
        # Get a list of the ensembles (could get this from the amopt dictionary)
        # For now we just use the submission scripts and assume all have .sh or .sub extension
        ext='.sh'
        if sys.platform.startswith("win"):
            ext='.bat'
        ensembles = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*"+ext) ) ]
        if not len(ensembles):
            # legacy - try .sub
            ensembles = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sub") ) ]

        if not len(ensembles):
            self.logger.warn("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            return False

        # reset any results
        self.results = []
        failed = {} # dict mapping failures to what went wrong - need to process at the end
        for ensemble in ensembles:

            # Check job directory
            jobDir = os.path.join( mrbumpDir, 'search_'+ensemble+'_mrbump' )
            #jobDir = os.path.join( mrbumpDir, 'search_'+ensemble )
            if not os.path.isdir(jobDir):
                self.logger.critical("Missing job directory: {0}".format( jobDir ) )
                failed[ ensemble ] = "no_job_directory"
                continue

            self.logger.debug(" -- checking directory for results: {0}".format( jobDir ) )
            # Check if finished
            if not os.path.exists( os.path.join( jobDir, "results", "finished.txt" ) ):
                self.logger.debug(" Found unfinished job: {0}".format( jobDir ) )
                failed[ ensemble ] = "unfinished"
                continue

            # Check resultsTable.dat
            resultsTable = os.path.join( jobDir,"results", "resultsTable.dat" )
            resultsDict = os.path.join( jobDir,"results", "resultsTable.pkl" )
            if os.path.isfile(resultsDict):
                self.mode="dict"
                self.results += self.processResultsPkl(resultsDict)
            elif os.path.isfile(resultsTable):
                self.results += self.parseTableDat(resultsTable)
            else:
                self.logger.debug(" -- Could not find results files: {0} or {1}".format(resultsTable,resultsDict))
                failed[ ensemble ] = "missing-results-file"
                continue

        # Process the failed results
        if failed and len(self.results):
            self._addFailedResults(mrbumpDir, failed)

        if not len(self.results):
            self.logger.warn("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            return False

        # Sort the results
        self.sortResultsD(self.results)

        return True

    def _getUnfinishedResult(self, result ):
        """Return a result for an unfinished job"""

        if not result.name:
            # Use directory name for job name
            dlist = os.listdir( os.path.join( result.jobDir, "data") )
            if len( dlist ) != 1:
                # something has gone really wrong...
                # Need to work out name from MRBUMP directory structure - search_poly_ala_trunc_6.344502_rad_3_phaser_mrbump
                dname = os.path.basename(result.jobDir)[7:-7]
                # Horrible - check if we were run with split_mr - in which case _phaser or _molrep are appended to the name
                if dname.endswith("_molrep") or dname.endswith("_phaser"):
                    dname = dname[:-7]
                # Add loc0_ALL_ and append _UNMOD. shudder...
                result.ensembleName = dname
                result.name = "loc0_ALL_" + dname + "_UNMOD"
                result.solution = "ERROR"
            else:
                # Use dirname but remove "loc0_ALL_" from front
                #result.name = os.path.basename( dlist[0] )[9:]
                # Use dirname but add "_UNMOD" to back
                result.ensembleName = os.path.basename( dlist[0] )
                result.name = result.ensembleName + "_UNMOD"

        result.program = "unknown"
        return

    def parseTableDat(self, tfile):
        
        tresults = self._parseTableDat(tfile)
        
        results=[]
        # Now convert to the dict form
        for r in tresults:
            d = {}
            d['Job_dir'] = r.jobDir
            d['MR_dir'] = r.mrDir
            d['Name'] = r.name
            d['SearchModel_name'] = r.ensembleName
            d['MR_program'] = r.program
            d['Solution_Type'] = r.solution
            d['final_Rfact'] = r.rfact
            d['final_Rfree'] = r.rfree
            d['SHELXE_CC'] = r.shelxeCC
            d['SHELXE_ACL'] = r.shelxeACL
            d['Bucc_final_Rfact'] = r.buccFinalRfact
            d['Bucc_final_Rfree'] = r.buccFinalRfree
            d['ARP_final_Rfact'] = r.arpWarpFinalRfact
            d['ARP_final_Rfree'] = r.arpWarpFinalRfree
            
            d['PHASER_pdb'] = r.phaserPdb
            d['MOLREP_pdb'] = r.molrepPdb
            d['REFMAC_pdb'] = r.refmacPdb
            d['SHELXE_pdb'] = r.shelxePdb
            d['Bucc_final_pdb'] = r.buccaneerPdb
            d['ARP_final_pdb'] = r.arpWarpPdb
            
            d['PHASER_TFZ'] = r.phaserTFZ
            d['PHASER_LLG'] = r.phaserLLG
            d['PHASER_time'] = r.phaserTime
            d['PHASER_killed'] = r.phaserKilled
            
            d['MOLREP_score'] = r.molrepScore
            d['MOLREP_time'] = r.molrepTime
            
            results.append(d)
        
        return results
    
    
    def _parseTableDat(self, tfile):
        """Read a resultsTable file and return a list of MrBump results objects"""

        # Extract the various components from the path
        tfile=os.path.abspath(tfile)
        paths = tfile.split( os.sep )
        jobDir = os.path.abspath( os.sep.join( paths[:-2] ) )
        ensemble = paths[-3][7:-7] #  remove search_..._mrbump: 'search_All_atom_trunc_0.005734_rad_1_mrbump'

        results = []
        header = None
        nfields=None
        # Read results table to get the results
        for line in open(tfile):

            # Create a result object for each line in the output file
            result = MrBumpResult()
            result.jobDir = jobDir
            result.ensembleName = ensemble

            line = line.strip()
            if not header:
                # Processing header
                header = line.split()
                nfields = len(header) # count as check
                for f in header:
                    # Map the data fields to their titles
                    if f not in self.title2attr.keys():
                        self.logger.critical("jobDir {0}: Problem with field {1} in headerline: {2}".format( jobDir, f, line ) )
                        result.header = ['Ensemble_Name']+[h for h in self.header if h in header ]# copy those attributes that we found
                        result.solution = "problem-header-file.dat"
                        self._getUnfinishedResult( result )
                        results.append( result )
                        return results
                continue
                # End header processing
            else:
                result.header = ['Ensemble_Name']+[h for h in self.header if h in header ] # copy those attributes that we found

            fields = line.split()
            if len(fields) != nfields:
                msg = "jobDir {0}: Problem getting dataline: {1}".format( jobDir, line )
                self.logger.debug(msg)
                result.solution = "corrupted-data-tfile.dat"
                self._getUnfinishedResult( result )
                results.append( result )
                continue

            # The headers tell us what attribute is in each column, so we use these and the header2attr dict to
            # set the results
            for i, title in enumerate(header):
                if title == 'MR_Program':
                    setattr( result, self.title2attr[title], fields[i].lower() )
                else:
                    setattr( result, self.title2attr[title], fields[i] )

            if result.program not in ['phaser','molrep']:
                raise RuntimeError,"getResult, unrecognised program in line: {0}".format(line)

            # Rebuild the path that generated the result
            # Add loc0_ALL_ and strip  _UNMOD from (e.g.): loc0_ALL_All_atom_trunc_0.34524_rad_1_UNMOD
            #dirName = "loc0_ALL_" + result.name + "_UNMOD"
            # While using old names - just strip _UNMOD
            dirName = result.name[:-6]
            result.mrDir = os.path.join( result.jobDir,'data',dirName,'unmod','mr',result.program )
            
            # See which pdb files were created
            if result.program=='phaser':
                phaserPdb=os.path.join(result.mrDir,
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result.program,result.ensembleName))
                if os.path.isfile(phaserPdb):
                    result.phaserPdb=phaserPdb
            elif result.program=='molrep':
                molrepPdb=os.path.join(result.mrDir,
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result.program, result.ensembleName) )
                if os.path.isfile(molrepPdb):
                    result.molrepPdb=molrepPdb

            refmacPdb=os.path.join(result.mrDir,
                                   'refine',
                                   "refmac_" + result.program + "_loc0_ALL_" + result.ensembleName + "_UNMOD.pdb" )
            if os.path.isfile(refmacPdb):
                result.refmacPdb=refmacPdb
                
            shelxePdb=os.path.join(result.mrDir,'build','shelxe',
                                   "shelxe_" + result.program + "_loc0_ALL_" + result.ensembleName + "_UNMOD.pdb" )
            if os.path.isfile(shelxePdb):
                result.shelxePdb=shelxePdb
                
            buccaneerPdb=os.path.join(result.mrDir,'build','shelxe','rebuild','buccaneer',
                                      "buccSX_output.pdb" )
            if os.path.isfile(buccaneerPdb):
                result.buccaneerPdb=buccaneerPdb
                
            arpWarpPdb=os.path.join(result.mrDir,'build','shelxe','rebuild','arpwarp',
                                      "refmacSX_output_warpNtrace.pdb" )
            if os.path.isfile(arpWarpPdb):
                result.arpWarpPdb=arpWarpPdb
            
            # Set the final result pdb
            if result.buccaneerPdb:
                result.pdb=result.buccaneerPdb
            elif result.arpWarpPdb:
                result.pdb=result.arpWarpPdb
            elif result.shelxePdb:
                result.pdb=result.shelxePdb
            elif result.refmacPdb:
                result.pdb=result.refmacPdb

            # Add the information from Shelxe
            #self.addShelxeResult( result )

            # Add the information from Buccaneer rebuild of the Shelxe trace
            #self.addBuccaneerResult( result )
            self.analyseResult(result)

            results.append( result )

        return results
    
    def processResultsPkl(self,resultsPkl):
        """Process dictionary of form:
        
        {'loc0_ALL_c1_tl100_r2_allatom_UNMOD':
           {'MOLREP': {'Molrep_time': None,
                       'Solution_Type': 'PHASER_FAIL',
                       'Phaser_killed': None,
                       'Phaser_time': None,
                       'ARP_final_Rfree': None,
                       'Molrep_score': None,
                       'Phaser_TFZ': None,
                       'ARP_final_Rfact': None,
                       'Bucc_final_Rfree': None,
                       'final_Rfree': 1.0,
                       'SHELXE_ACL': None,
                       'final_Rfact': 1.0,
                       'SHELXE_CC': None,
                       'SHELXE_time': None,
                       'SHELXE_MCL': None,
                       'Bucc_final_Rfact': None,
                       'SHELXE_NC': None,
                       'Phaser_LLG': None},
            {'PHASER' : ...
        """
        with open(resultsPkl) as f:
            rD=cPickle.load(f)
        if not rD: return []
        
        results=[]
        for name,d1 in rD.iteritems():
            for mrprog,d2 in d1.iteritems():
                # Add MR program as dictionary entry
                d = copy.copy(d2)
                d['SearchModel_name'] = name
                d['MR_program'] = mrprog
                results.append(d)
                
        return results
    
    def analyseResult(self,result):
        
        mrDir = result.mrDir
        
        result.ensembleName = result.name[9:-6]
        
        # Process log
        #mrbumpP = MrbumpLogParser(result.mrbumpLog)
        #result.estChainsASU=mrbumpP.noChainsTarget
        
        if result.program == "phaser":
            if result.phaserPdb:
                phaserP = parse_phaser.PhaserPdbParser(result.phaserPdb)
                result.phaserLLG = phaserP.LLG
                result.phaserTFZ = phaserP.TFZ
                result.phaserPdb = result.phaserPdb
            
            phaserLog = os.path.join( mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(result.program, result.ensembleName) )
            if os.path.isfile( phaserLog ):
                phaserP = parse_phaser.PhaserLogParser( phaserLog, onlyTime=True )
                result.phaserLog    = phaserLog
                result.phaserTime   = phaserP.time
                result.phaserKilled = phaserP.killed
            
        elif result.program == "molrep":
            molrepLog = os.path.join( mrDir, "molrep.log" )
            result.molrepLog = molrepLog
            molrepP = MolrepLogParser( molrepLog )
            result.molrepScore = molrepP.score
            result.molrepTime = molrepP.time
        #
        # SHELXE PROCESSING
        #
# We get this from the pdb
#         shelxeLog = os.path.join( mrDir, "build/shelxe/shelxe_run.log" )
#         if os.path.isfile( shelxeLog ):
#             result.shelxeLog = shelxeLog
#             shelxeP = parse_shelxe.ShelxeLogParser( shelxeLog )
#             result.shelxeCC = shelxeP.CC
#             result.shelxeACL = shelxeP.avgChainLength
#             result.shelxeMaxChainLength = shelxeP.maxChainLength
#             result.shelxeNumChains= shelxeP.numChains
        #
        # Buccaneer Rebuild Processing
        #
        buccaneerLog = os.path.join( mrDir,
                                     "build/shelxe/rebuild/buccaneer",
                                     "buccaneer.log" )
    
        bp = parse_buccaneer.BuccaneerLogParser()
        if os.path.isfile( buccaneerLog ):
            bp.parse( buccaneerLog )
            result.buccFinalRfree = bp.finalRfree
            result.buccFinalRfact = bp.finalRfact

        #
        # Arpwarp Rebuild Processing
        #
        arpLog = os.path.join(mrDir,
                              "build/shelxe/rebuild/arpwarp",
                              "arpwarp.log")
        if os.path.isfile(arpLog):
                        ap=parse_arpwarp.ArpwarpLogParser()
                        ap.parse(arpLog)
                        result.arpWarpFinalRfact=ap.finalRfact
                        result.arpWarpFinalRfree=ap.finalRfree
        
        return

    def sortResults( self, results ):
        """
        Sort the results
        """
        
        # Check each result to see what attributes are set and use this to work out how we rate this run
        reverse=False
        sortf=False
        for r in results:
            if r.shelxeCC and r.shelxeCC != "--" and float(r.shelxeCC) > -1.0:
                reverse=True
                sortf = lambda x: float( x.shelxeCC )
                break
            elif r.buccRfree and r.buccRfree != "--" and float(r.buccRfree) < 1.0:
                sortf = lambda x: float( x.buccRfree )
                break
            elif r.arpWarpRfree and r.arpWarpRfree != "--" and float(r.arpWarpRfree) < 1.0:
                sortf = lambda x: float( x.arpWarpRfree )
                break
            elif r.rfree and float(r.rfree) < 1.0:
                sortf = lambda x: float( x.rfree )
                break

        if sortf:
            # Now sort by the key
            results.sort(key=sortf, reverse=reverse)
        return
    
    def sortResultsD( self, results ):
        """
        Sort the results
        """
        
        # Check each result to see what attributes are set and use this to work out how we rate this run
        reverse=False
        sortf=False
        for r in results:
            if 'SHELXE_CC' in r and r['SHELXE_CC'] and float(r['SHELXE_CC']) > -1.0:
                reverse=True
                sortf = lambda x: float(0) if x['SHELXE_CC']  is None else float( x['SHELXE_CC'])
                break
            # Not sure what to do about buccaneer/arpwarp when no shelxe?
#             elif r.buccRfree and r.buccRfree != "--" and float(r.buccRfree) < 1.0:
#                 sortf = lambda x: float( x.buccRfree )
#                 break
#             elif r.arpWarpRfree and r.arpWarpRfree != "--" and float(r.arpWarpRfree) < 1.0:
#                 sortf = lambda x: float( x.arpWarpRfree )
#                 break
            elif 'final_Rfree' in r and r['final_Rfree'] and float(r['final_Rfree']) < 1.0:
                sortf = lambda x: float('inf') if x['final_Rfree']  is None else float( x['final_Rfree'] )
                break

        if sortf:
            # Now sort by the key
            results.sort(key=sortf, reverse=reverse)
        return

    def summariseResults( self, mrbumpDir ):
        """Return a string summarising the results"""

        got = self.extractResults( mrbumpDir )
        if got:
            return self.summaryString()
        else:
            return "\n!!! No results found in directory: {0}\n".format( mrbumpDir )

    def summaryString( self ):
        """Return a string suitable for printing the sorted results"""

        resultsTable = []
        keys = ['SearchModel_name','MR_program','Solution_Type','final_Rfact','final_Rfree','SHELXE_CC',
                'SHELXE_ACL','Bucc_final_Rfact','Bucc_final_Rfree','ARP_final_Rfact','ARP_final_Rfree']
        resultsTable.append(keys)
        for result in self.results:
            resultLine = []
            for k in keys:
                resultLine.append(result[k])
            resultsTable.append( resultLine )

        # Format the results
        table = printTable.Table()
        summary = table.pprint_table( resultsTable )

        r = "\n\nOverall Summary:\n\n"
        r += summary

#         r += '\nBest Molecular Replacement results so far are in:\n\n'
#         r +=  self.results[0]['MR_dir']

#         if self.results[0].pdb and os.path.isfile(self.results[0].pdb):
#             r += '\n\nFinal PDB is:\n\n'
#             r +=  self.results[0].pdb
        r += '\n\n'

        return r

class Test(unittest.TestCase):


    def testResultsDict(self):
        """Parse a results file"""
        
        resultsPkl="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_0/MRBUMP/MRBUMP/search_c1_tl100_r2_allatom_mrbump/results/resultsTable.pkl"

        rs=ResultsSummary()
        print rs.processResultsPkl(resultsPkl)
        
    def testProcess(self):
        """Parse a results file"""
        
        mrbdir="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_0/MRBUMP/MRBUMP"
        rs=ResultsSummary()
        print rs.summariseResults(mrbdir)

        return

if __name__ == "__main__":
    if len(sys.argv) == 2:
        mrbump_dir = os.path.join( os.getcwd(), sys.argv[1] )
    else:
        mrbump_dir = "/Users/jmht/Documents/AMPLE/res.test/cluster_run1"
        mrbump_dir = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/MRBUMP/cluster_1"
        mrbump_dir = "/gpfs/home/HCEA041/djr01/jxt15-djr01/TM/3OUF/ROSETTA_MR_1/MRBUMP/cluster_run1"

    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    r = ResultsSummary()
    print r.summariseResults( mrbump_dir )

