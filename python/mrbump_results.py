#!/usr/bin/env ccp4-python

import copy
import glob
import logging
import os
import types

# Hack to make sure we can find the modules we need
if __name__ == "__main__":
    import sys
    root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
    sys.path.insert( 0, os.path.join( root, "scripts" ) )


# Our imports
import parse_buccaneer
import parse_shelxe
import printTable

# We need a null logger so that we can be used without requiring a logger
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

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
        self.shelxCC = -1.0
        self.shelxeAvgChainLength = -1.0
        # Buccaneer scores after final rebuild
        self.buccFinalRfact = 1.0
        self.buccFinalRfree  = 1.0
        
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
                            'Bucc_final_Rfact' : 'buccRfact',
                            'Bucc_final_Rfree' : 'buccRfree',
                            'ARP_final_Rfact'  : 'arpWarpRfact',
                            'ARP_final_Rfree'  : 'arpWarpRfree',
                            'SHELXE_CC'        : 'shelxCC',
                            'SHELXE_ACL'       : 'shelxeAvgChainLength',
                             }
        
        self.logger = logging.getLogger()
        # Add Null logger so we can be used without requiring a logger
        self.logger.addHandler(NullHandler())
        return
        
    def _addFailedResults(self, mrbumpDir, failed, header):
        """Add failures to self.results
        
        Args:
        failed: dict of {ensemble : result}
        header: list with header for results table
        """
        assert failed and header
        
        count=0
        for ensemble, reason in failed.iteritems():
            result = MrBumpResult()
            # name hard-coded
            result.name = "loc0_ALL_" + ensemble + "_UNMOD"
            result.jobDir = os.path.join( mrbumpDir, 'search_'+ensemble+'_mrbump' )
            result.header = header
            result.solution = reason
            self._getUnfinishedResult( result )
            self.results.append( result )
            count += 1
        
        self.logger.debug("Added {0} MRBUMP result failures".format(count) )
        return
    
    def addBuccaneerResult(self, result ):
        """Add the buccaneer rebuild information to the result"""


        buccaneerPdb = os.path.join( result.mrDir,
                                     "build/shelxe/rebuild/build",
                                     "buccSX_output.pdb" )
        buccaneerLog = os.path.join( result.mrDir,
                                     "build/shelxe/rebuild/build",
                                     "buccaneer.log" )
        
        bp = parse_buccaneer.BuccaneerLogParser()
        if os.path.isfile( buccaneerLog ):
            bp.parse( buccaneerLog )
            result.buccFinalRfree = bp.finalRfree
            result.buccFinalRfact = bp.finalRfact
        
        if os.path.isfile( buccaneerPdb ):
            result.buccaneerPdb = buccaneerPdb

        return

    def addShelxeResult(self, result ):
        """Add the shelxe information to the result"""
        
        shelxePdb = os.path.join( result.mrDir,
                                  "build",
                                  "shelxe",
                                  "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( result.program,
                                                                              result.ensembleName ) )
        if os.path.isfile( shelxePdb):
            result.shelxePdb = shelxePdb
            
        shelxeLog = os.path.join( result.mrDir,
                                  "build",
                                  "shelxe",
                                  "shelxe_run.log" )
        
        if os.path.isfile( shelxeLog ):
            shelxeP = parse_shelxe.ShelxeLogParser( shelxeLog )
            #assert result.shelxCC == shelxeP.CC,"Mismatching ShelxeCC scores"
            result.shelxeAvgChainLength = shelxeP.avgChainLength
            result.shelxeLog            = shelxeLog
            result.shelxeMaxChainLength = shelxeP.maxChainLength
            result.shelxeNumChains      = shelxeP.numChains

            # Another horrible hack - add the title to the header
            result.header.append('SHELXE_ACL')
        else:
            #assert result.shelxCC == shelxeP.CC,"Mismatching ShelxeCC scores"
            result.shelxeAvgChainLength = 0
            result.shelxeLog            = ""
            result.shelxeMaxChainLength = 0
            result.shelxeNumChains      = 0

            # Another horrible hack - add the title to the header
            result.header.append('SHELXE_ACL')

        return


    def extractResults( self, mrbumpDir ):
        """
        Find the results from running MRBUMP and sort them
        """
        mrbumpDir = os.path.abspath( mrbumpDir )
        # Get a list of the ensembles (could get this from the amopt dictionary)
        # For now we just use the submission scripts and assume all have .sh or .sub extension
        ensembles = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sh") ) ]
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
            if not os.path.exists(resultsTable):
                self.logger.debug(" -- Could not find file: {0}".format( resultsTable ) )
                failed[ ensemble ] = "missing-resultsTable.dat"
                continue
            
            # Extract the result
            self.results += self.parseTableDat(resultsTable)
            
#         if not header or not len(header):
#             self.logger.warn("Could not extract any results from directory - no header: {0}".format( mrbumpDir ) )
#             return False

        # Process the failed results
        if failed and len(self.results):
            self._addFailedResults( mrbumpDir, failed, self.results[0].header )
                
        if not len(self.results):
            self.logger.warn("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            return False

        # Sort the results
        self.sortResults(self.results)
        
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
        """Read a resultsTable file and return a list of MrBump results objects"""
        
        # Extract the various components from the path
        paths = tfile.split( os.sep )
        assert len( paths[0] )  == 0,"Need absolute path: {0}".format( tfile )
        
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
                        result.header = header
                        result.solution = "problem-header-file.dat"
                        self._getUnfinishedResult( result )
                        results.append( result )
                        return results
                continue
                # End header processing
            else:
                # horrible - we manipulate the header in addShelxe/buccaneer so we need to use a copy here
                result.header = copy.copy(header)

            
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
            
            # Need to reconstruct final pdb from directory and program name
            pdbName = "refmac_" + result.program + "_loc0_ALL_" + result.ensembleName + "_UNMOD.pdb"
            result.pdb = os.path.join( result.mrDir,'refine', pdbName )
            
            # Add the information from Shelxe
            #self.addShelxeResult( result )
            
            # Add the information from Buccaneer rebuild of the Shelxe trace
            #self.addBuccaneerResult( result )
            
            results.append( result )
            
        return results
        
    def sortResults( self, results ):
        """
        Sort the results
        """
        
        # Check each result to see what attributes are set and use this to work out how we rate this run
        reverse=False
        sortf=False
        for r in results:
            if r.shelxCC and r.shelxCC != "--" and float(r.shelxCC) > -1.0:
                reverse=True
                sortf = lambda x: float( x.shelxCC )
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
        
        #Header
        # The best result may have more info (e.g. shelxe info) aso we use this to 
        # work out what data we should present
        topHeader = self.results[0].header
        resultsTable.append( topHeader )
        
        for result in self.results:
            resultLine = []
            #for h in result.header:
            for h in topHeader:
                resultLine.append( getattr( result, self.title2attr[h] ) )
            resultsTable.append( resultLine )

        # Format the results
        table = printTable.Table()
        summary = table.pprint_table( resultsTable )
        
        r = "\n\nOverall Summary:\n\n"
        r += summary
        r += '\nBest results so far are in :\n\n'
        r +=  self.results[0].mrDir
            
        return r    
#
# DEPRECATED CODE BELOW HERE
#
#def get_rfree(pdb):
#  pdb = open(pdb)
#
#  freer_pattern = re.compile('REMARK\s*\d*\s*FREE R VALUE\s*\:\s*(\d*\.\d*)')
#  FreeR = 'nan' 
#
#  for line in pdb:
#
#    if re.search(freer_pattern, line):
#    #  print line
#      split=  re.split(freer_pattern, line)
#   #   print split
#      FreeR = split[1]
#
#  return FreeR
############################
#
#def get_shelx_score(RESULT):
#
#    Best_result = 0
#    result = open(RESULT)
#    fail_pattern = ('\*\* Unable to trace map - giving up \*\*')
#    pattern = re.compile('CC for partial structure against native data =\s*(\d*\.\d*)\s*%')
#    for line in result: 
#       if re.search(pattern, line):
#        split= re.split(pattern, line)
#      #  print split
#        if float(split[1]) > Best_result:
#          Best_result = float(split[1])
#       if re.search(fail_pattern, line):
#        Best_result = line.rstrip('\n')
#
#    return Best_result
#  
###########################
#def make_log(mr_bump_path, final_log):
#  final_log = open(final_log, "w")
#
#  list_dir = os.listdir(mr_bump_path)
#  for a_dir in list_dir:
#    if os.path.isdir(mr_bump_path + '/'+a_dir):
#
#     name=re.sub('search_', '', a_dir)
#     name=re.sub('_mrbump', '', name)
#     
#     phaser_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'
#     molrep_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'
#     phaser_log = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/phaser/phaser_loc0_ALL_'+name+'_UNMOD.1.pdb'
#
#     shelx_phaser = mr_bump_path + '/'+a_dir+'/phaser_shelx/RESULT'
#     shelx_molrep = mr_bump_path + '/'+a_dir+'/molrep_shelx/RESULT'
#
#     if os.path.exists(phaser_pdb):
#      if os.path.exists(shelx_phaser):
#        phaser_FreeR =   get_rfree(phaser_pdb)  
#        phaser_shelx =   get_shelx_score(shelx_phaser)
#        final_log.write('Ensembe ' + name+'  phaser FreeR: '+phaser_FreeR +  '  shelx score: '+str(phaser_shelx) + '\n')
#        final_log.flush()
#  
#     if os.path.exists(molrep_pdb):
#      if os.path.exists(shelx_molrep):
#        molrep_FreeR =   get_rfree(molrep_pdb) 
#        molrep_shelx =   get_shelx_score(shelx_molrep)  
#        final_log.write('Ensembe ' + name+'  molrep FreeR: '+molrep_FreeR +  '  shelx score: '+str(molrep_shelx) + '\n')
#        final_log.flush()
#
#
#def get_shel_score(pdb):
#  print 'here'
#  pdb=open(pdb)
#  score = 0
#
#  pattern=re.compile('CC\s*=\s*(\d*\.\d*)')
#  for line in pdb:
#  
#   if re.search(pattern, line):
#     
#     split = re.split(pattern,line)
#     
#     score = split[1]
#  return score
## # # # # # # # # # # # # # # # 
#
#def rank_results(mrbump_path, overpath):
# 
# order = []
#
# if not os.path.exists( os.path.join(overpath, 'RESULTS')): 
#   os.mkdir(os.path.join(overpath, 'RESULTS'))
#
# 
#
# for folder in os.listdir(mrbump_path):
#  if  os.path.isdir(mrbump_path+'/'+folder):
#    
#
#    if os.path.exists(mrbump_path+'/'+folder+'/phaser_shelx'): # phaser shelx
#       if os.path.exists(mrbump_path+'/'+folder+'/phaser_shelx/orig.pdb'):
#          score = get_shel_score(mrbump_path+'/'+folder+'/phaser_shelx/orig.pdb')       
#          order.append([mrbump_path+'/'+folder+'/phaser_shelx/XYZOUT', float( score), mrbump_path+'/'+folder+'/phaser_shelx/HKLOUT' ])
#
#    if os.path.exists(mrbump_path+'/'+folder+'/molrep_shelx'): # phaser shelx
#       if os.path.exists(mrbump_path+'/'+folder+'/molrep_shelx/orig.pdb'):
#          score = get_shel_score(mrbump_path+'/'+folder+'/molrep_shelx/orig.pdb')
#          order.append([mrbump_path+'/'+folder+'/molrep_shelx/XYZOUT', float( score), mrbump_path+'/'+folder+'/molrep_shelx/HKOUT' ])
#          print 'GOT'
#
# 
#
# index=1
# order.sort(key=lambda tup: tup[1])
# order.reverse()
# for x in order:
#   print x
#   os.system('cp  ' + x[0] + ' '+ os.path.join(overpath, 'RESULTS', 'result_'+str(index)+'.pdb') )
#   os.system('cp  ' + x[2] + ' '+ os.path.join(overpath, 'RESULTS', 'result_'+str(index)+'.mtz') )
# 
#   index+=1
#############
##mr_bump_path = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP'
##final_log = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP/FINAL'
##make_log(mr_bump_path,final_log)


if __name__ == "__main__":
    import sys
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
    
