#!/usr/bin/env ccp4-python

import copy
import cPickle
import glob
import logging
import os
import sys
import unittest

# Hack to make sure we can find the modules we need
if __name__ == "__main__":
    root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
    sys.path.insert( 0, os.path.join( root, "scripts" ) )

# Our imports
if not "CCP4" in os.environ.keys():
    raise RuntimeError('CCP4 not found')
mrbumpd=os.path.join(os.environ['CCP4'],"include","mrbump","include","parsers")
mrbumpd="/opt/mrbump-trunk/include/parsers"
sys.path.insert(0,mrbumpd)
import parse_arpwarp
import parse_buccaneer
import parse_phaser
import printTable

# We need a null logger so that we can be used without requiring a logger
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """

    def __init__(self):

        self.results = []
        # List of all the possible column titles and their result object attributes
        self.title2key = {
                            'Model_Name'       : 'name',
                            'MR_Program'       : 'MR_program',
                            'Solution_Type'    : 'Solution_Type',
                            'final_Rfact'      : 'final_Rfact',
                            'final_Rfree'      : 'final_Rfree',
                            'SHELXE_CC'        : 'SHELXE_CC',
                            'SHELXE_ACL'       : 'SHELXE_ACL',
                            'Bucc_final_Rfact' : 'BUCC_final_Rfact',
                            'Bucc_final_Rfree' : 'BUCC_final_Rfree',
                            'ARP_final_Rfact'  : 'ARP_final_Rfact',
                            'ARP_final_Rfree'  : 'ARP_final_Rfree',
                             }

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
            #d['name'] = "loc0_ALL_" + ensemble + "_UNMOD"
            d['name'] = "loc0_ALL_" + ensemble + "_UNMOD"
            d['ensemble_name']=ensemble 
            d['JobDirectory'] = os.path.join( mrbumpDir, 'search_'+ensemble+'_mrbump' )
            d['Solution_Type'] = reason
            self.results.append( d )
            count += 1
        self.logger.debug("Added {0} MRBUMP result failures".format(count) )
        return

    def createDict(self):
        d = {}
        
        # our additional keys
        d['ensemble_name'] = None
        d['MR_program'] = None
        d['name'] = None
        
        d['JobDirectory']        = None
        d['Solution_Type']       = None
        
        d['PHASER_LLG']         = None
        d['PHASER_TFZ']         = None
        d['PHASER_RFZ']         = None
        d['PHASER_time']        = None
        d['PHASER_killed']      = None
        d['PHASER_pdbout']      = None
        d['PHASER_mtzout']      = None
        d['PHASER_logfile']      = None
        d['PHASER_version']         = None
        
        d['MOLREP_score']      = None
        d['MOLREP_time']      = None
        d['MOLREP_pdbout']      = None
        d['MOLREP_logfile']      = None
        d['MOLREP_version']         = None
        
        d['final_Rfact']       = None
        d['final_Rfree']       = None
        d['REFMAC_pdbout']      = None
        d['REFMAC_mtzout']      = None
        d['REFMAC_logfile']      = None
        d['REFMAC_version']         = None
        
        d['BUCC_final_Rfact']       = None
        d['BUCC_final_Rfree']       = None
        d['BUCC_pdbout']      = None
        d['BUCC_mtzout']      = None
        d['BUCC_logfile']      = None
        d['BUCC_version']         = None
        
        d['ARP_final_Rfact']       = None
        d['ARP_final_Rfree']       = None
        d['ARP_pdbout']      = None
        d['ARP_mtzout']      = None
        d['ARP_logfile']      = None
        d['ARP_version']         = None
        
        d['SHELXE_CC']          = None
        d['SHELXE_ACL']         = None
        d['SHELXE_MCL']         = None
        d['SHELXE_NC']         = None
        d['SHELXE_wMPE']     = None
        d['SHELXE_os']         = None
        d['SHELXE_time']     = None
        d['SHELXE_pdbout']     = None
        d['SHELXE_phsout']     = None
        d['SHELXE_mtzout']     = None
        d['SHELXE_logfile']     = None
        d['SHELXE_version']     = None
        
        d['SXRBUCC_version']        = None
        d['SXRBUCC_final_Rfact']    = None
        d['SXRBUCC_final_Rfree']    = None
        d['SXRBUCC_pdbout']      = None
        d['SXRBUCC_mtzout']      = None
        d['SXRBUCC_logfile']      = None
        
        d['SXRARP_version']         = None
        d['SXRARP_final_Rfact']     = None
        d['SXRARP_final_Rfree']     = None
        d['SXRARP_pdbout']      = None
        d['SXRARP_mtzout']      = None
        d['SXRARP_logfile']      = None
        
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
            resultsDict = os.path.join( jobDir,"results", "resultsTable.pkl" )
            resultsTable = os.path.join( jobDir,"results", "resultsTable.dat" )
            if os.path.isfile(resultsDict):
                self.results += self.processResultsPkl(resultsDict)
            elif os.path.isfile(resultsTable):
                self.results += self.parseTableDat(resultsTable)
            else:
                self.logger.debug(" -- Could not find results files: {0} or {1}".format(resultsDict,resultsTable))
                failed[ ensemble ] = "missing-results-file"
                continue

        # Process the failed results
        if failed and len(self.results):
            self._addFailedResults(mrbumpDir, failed)

        if not len(self.results):
            self.logger.warn("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            return False

        # Sort the results
        self.sortResults(self.results)
        return True

    def _getUnfinishedResult(self, result ):
        """Return a result for an unfinished job"""
        if not result['name']:
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
                result['ensemble_name'] = dname
                result['name'] = "loc0_ALL_" + dname + "_UNMOD"
                result["Solution_Type"] = "ERROR"
            else:
                # Use dirname but remove "loc0_ALL_" from front
                #result.name = os.path.basename( dlist[0] )[9:]
                # Use dirname but add "_UNMOD" to back
                result['ensemble_name'] = os.path.basename( dlist[0] )
                result['name'] = result.ensembleName + "_UNMOD"
        result["MR_program"] = "unknown"
        return
    
    def parseTableDat(self, tfile):
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
            result = self.createDict()
            result['ensemble_name'] = ensemble

            line = line.strip()
            if not header:
                # Processing header
                header = line.split()
                nfields = len(header) # count as check
                for f in header:
                    # Map the data fields to their titles
                    if f not in self.title2key.keys():
                        self.logger.critical("jobDir {0}: Problem with field {1} in headerline: {2}".format( jobDir, f, line ) )
                        result['Solution_Type'] = "problem-header-file.dat"
                        self._getUnfinishedResult( result )
                        results.append( result )
                        return results
                continue
                # End header processing

            fields = line.split()
            if len(fields) != nfields:
                msg = "jobDir {0}: Problem getting dataline: {1}".format( jobDir, line )
                self.logger.debug(msg)
                result['Solution_Type'] = "corrupted-data-tfile.dat"
                self._getUnfinishedResult( result )
                results.append( result )
                continue

            # The headers tell us what attribute is in each column, so we use these and the header2attr dict to
            # set the results
            for i, title in enumerate(header):
                v=fields[i]
                if v == '--': v=None # non-valid values in table are indicated by --
                result[self.title2key[title]] = v

            dirName = result['name'][:-6]
            result["JobDirectory"] = os.path.join( jobDir,'data',dirName,'unmod','mr',result["MR_program"].lower() )
            
            # See which pdb files were created
            if result["MR_program"]=='PHASER':
                phaserPdb=os.path.join(result["JobDirectory"],
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result["MR_program"].lower(),result['ensemble_name']))
                if os.path.isfile(phaserPdb):
                    result.phaserPdb=phaserPdb
            elif result.program=='molrep':
                molrepPdb=os.path.join(result["JobDirectory"],
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result["MR_program"].lower(),result['ensemble_name']) )
                if os.path.isfile(molrepPdb):
                    result.molrepPdb=molrepPdb

            refmacPdb=os.path.join(result["JobDirectory"],
                                   'refine',
                                   "refmac_" + result["MR_program"].lower() + "_loc0_ALL_" + result['ensemble_name'] + "_UNMOD.pdb" )
            if os.path.isfile(refmacPdb):
                result["REFMAC_pdbout"]=refmacPdb
                
            shelxePdb=os.path.join(result["JobDirectory"],'build','shelxe',
                                   "shelxe_" + result["MR_program"].lower() + "_loc0_ALL_" + result['ensemble_name'] + "_UNMOD.pdb" )
            if os.path.isfile(shelxePdb):
                result["SHELXE_pdbout"]=shelxePdb
                
            buccaneerPdb=os.path.join(result["JobDirectory"],'build','shelxe','rebuild','buccaneer',
                                      "buccSX_output.pdb" )
            if os.path.isfile(buccaneerPdb):
                result["SXRBUCC_pdbout"]=buccaneerPdb
                
            arpWarpPdb=os.path.join(result["JobDirectory"],'build','shelxe','rebuild','arpwarp',
                                      "refmacSX_output_warpNtrace.pdb" )
            if os.path.isfile(arpWarpPdb):
                result["SXRARP_pdbout"]=arpWarpPdb
            
            self.analyseResult(result)

            results.append( result )

        return results

    def processResultsPkl(self,resultsPkl):
        """Process dictionary
        """
        with open(resultsPkl) as f:
            rD=cPickle.load(f)
        if not rD: return []
        
        results=[]
        for name,d1 in rD.iteritems():
            for mrprog,d2 in d1.iteritems():
                # Check if all the entries are None - means job didn't run.
                # Should probably think of a better way to spot that
                if not any(d2.values()): continue
                # Add MR program as dictionary entry
                d = copy.copy(d2)
                del d['SearchModel_filename']
                d['name'] = name
                # name is e.g.: loc0_ALL_c1_tl100_r2_allatom_UNMOD
                d['ensemble_name'] = name[9:-6]
                d['MR_program'] = mrprog
                results.append(d)
        return results
    
    def analyseResult(self,result):
        
        mrDir = result["JobDirectory"]
        
        #result.ensembleName = result.name[9:-6]
        
        # Process log
        #mrbumpP = MrbumpLogParser(result.mrbumpLog)
        #result.estChainsASU=mrbumpP.noChainsTarget
        
        if result["MR_program"] == "PHASER":
            if result["PHASER_pdbout"]:
                phaserP = parse_phaser.PhaserPdbParser(result["PHASER_pdbout"])
                result["PHASER_LLG"] = phaserP.LLG
                result["PHASER_TFZ"] = phaserP.TFZ
            
            phaserLog = os.path.join( mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(result["MR_program"].lower(),result['ensemble_name']) )
            if os.path.isfile( phaserLog ):
                phaserP = parse_phaser.PhaserLogParser( phaserLog, onlyTime=True )
                #result.phaserLog    = phaserLog
                result["PHASER_time"]   = phaserP.time
                result["PHASER_killed"] = phaserP.killed
            
        #
        # SHELXE PROCESSING
        #
        #
        # Buccaneer Rebuild Processing
        #
        buccaneerLog = os.path.join( mrDir,
                                     "build/shelxe/rebuild/buccaneer",
                                     "buccaneer.log" )
    
        bp = parse_buccaneer.BuccaneerLogParser()
        if os.path.isfile( buccaneerLog ):
            bp.parse( buccaneerLog )
            result["SXRBUCC_final_Rfree"] = bp.finalRfree
            result["SXRBUCC_final_Rfact"] = bp.finalRfact

        #
        # Arpwarp Rebuild Processing
        #
        arpLog = os.path.join(mrDir,
                              "build/shelxe/rebuild/arpwarp",
                              "arpwarp.log")
        if os.path.isfile(arpLog):
                        ap=parse_arpwarp.ArpwarpLogParser()
                        ap.parse(arpLog)
                        result["SXRARP_final_Rfact"]=ap.finalRfact
                        result["SXRARP_final_Rfree"]=ap.finalRfree
        
        return

    def sortResults( self, results ):
        """
        Sort the results
        """
        # Check each result to see what attributes are set and use this to work out how we rate this run
        reverse=False
        sortf=False
        for r in results:
            if 'SHELXE_CC' in r and r['SHELXE_CC'] and float(r['SHELXE_CC']) > 0.0:
                reverse=True
                sortf = lambda x: float(0) if x['SHELXE_CC']  is None else float( x['SHELXE_CC'])
                break
            elif 'BUCC_final_Rfact' in r and r['BUCC_final_Rfact'] and float(r['BUCC_final_Rfact']) < 1.0:
                sortf = lambda x: float('inf') if x['BUCC_final_Rfact']  is None else float( x['BUCC_final_Rfact'] )
                break
            elif 'ARP_final_Rfree' in r and r['ARP_final_Rfree'] and float(r['ARP_final_Rfree']) < 1.0:
                sortf = lambda x: float('inf') if x['ARP_final_Rfree']  is None else float( x['ARP_final_Rfree'] )
                break
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
        keys = ['ensemble_name','MR_program','Solution_Type','final_Rfact','final_Rfree','SHELXE_CC',
                'SHELXE_ACL','BUCC_final_Rfact','BUCC_final_Rfree','ARP_final_Rfact','ARP_final_Rfree']
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

        r += '\nBest Molecular Replacement results so far are in:\n\n'
        r +=  self.results[0]["JobDirectory"]

#         if self.results[0].pdb and os.path.isfile(self.results[0].pdb):
#             r += '\n\nFinal PDB is:\n\n'
#             r +=  self.results[0].pdb
        r += '\n\n'

        return r
    
def checkSuccess(script_path):
    """
    Check if a job ran successfully.
    
    Args:
    directory -- directory mr bump ran the job
    
    Returns:
    True if success
    
    Success is assumed as a SHELX CC score of >= SHELXSUCCESS
    """
    directory, script = os.path.split(script_path)
    scriptname = os.path.splitext( script )[0]
    rfile = os.path.join(directory, 'search_'+scriptname+'_mrbump','results/resultsTable.pkl')
    #print "{0} checking for file: {1}".format(multiprocessing.current_process().name,rfile)
    if not os.path.isfile(rfile):
        #print "{0} cannot find results file: {1}".format(multiprocessing.current_process().name,rfile)
        return False
    
    # Results summary object to parse table file
    mrbR = ResultsSummary()
    
    # Put into order and take top one
    results = mrbR.processResultsPkl(rfile)
    mrbR.sortResults(results)
    r = results[0]

    success=False
    rFreeSuccess=0.4
    if 'SHELXE_CC' in r and r['SHELXE_CC'] and float(r['SHELXE_CC']) >= 25.0:
        success=True
    elif 'BUCC_final_Rfact' in r and r['BUCC_final_Rfact'] and float(r['BUCC_final_Rfact']) <= rFreeSuccess:
        success=True
    elif 'ARP_final_Rfree' in r and r['ARP_final_Rfree'] and float(r['ARP_final_Rfree']) <= rFreeSuccess:
        success=True
    elif 'final_Rfree' in r and r['final_Rfree'] and float(r['final_Rfree']) <= rFreeSuccess:
        success=True
        
    return success
    

def finalSummary(amoptd):
    """Print a final summary of the job"""
    
    
    ensembles_data=amoptd['ensembles_data']
    mrbump_data=amoptd['mrbump_results']
    
    # Merge dictionaries together
    results=[]
    for mrb in mrbump_data:
        d=copy.copy(mrb)
        for ed in ensembles_data:
            if ed['name'] == d['ensemble_name']:
                d.update(ed)
                results.append(d)

    resultsTable = []
    keys = ['ensemble_name','MR_program',"PHASER_LLG","PHASER_TFZ",'SHELXE_CC','SHELXE_ACL',"SXRBUCC_final_Rfact","SXRBUCC_final_Rfree",
            'SXRARP_final_Rfact','SXRARP_final_Rfree','subcluster_num_models','truncation_num_residues']

    resultsTable.append(keys)
    for result in results:
        resultLine = []
        for k in keys:
            resultLine.append(result[k])
        resultsTable.append( resultLine )

    # Format the results
    table = printTable.Table()
    summary = table.pprint_table( resultsTable )

    r = "\n\nOverall Summary:\n\n"
    r += summary

    r += '\nBest Molecular Replacement results so far are in:\n\n'
    r +=  results[0]["JobDirectory"]
    r += '\n\n'
    return r


class Test(unittest.TestCase):

    def testResultsDict(self):
        """Parse a results file"""
        
        resultsPkl="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_4/MRBUMP/MRBUMP/search_c1_tl100_r2_allatom_mrbump/results/resultsTable.pkl"

        rs=ResultsSummary()
        print rs.processResultsPkl(resultsPkl)
        
    def testProcess(self):
        """Parse a results file"""
        
        mrbdir="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_0/MRBUMP/MRBUMP"
        rs=ResultsSummary()
        print rs.summariseResults(mrbdir)
        return
    
    def testFinalSummary(self):
        """Parse a results file"""
        pkl="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_4/resultsd.pkl"
        with open(pkl) as f:
            d=cPickle.load(f)
        print finalSummary(d)
        return

if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print "Usage: {0} <MRBUMP_directory>".format(sys.argv[0])
        
    mrbump_dir = os.path.join( os.getcwd(), sys.argv[1] )
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    r = ResultsSummary()
    print r.summariseResults( mrbump_dir )

