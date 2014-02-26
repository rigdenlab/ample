#!/usr/bin/env python

import glob
import locale
import logging
import os
import re
import types

# Our imports
import printTable

class MrBumpResult(object):
    """
    Class to hold the result of running a MRBUMP job
    """
    def __init__(self):
        """
        
        """
        self.jobDir = None # directory jobs ran in
        self.resultDir = None # where the actual results are
        self.pdb = None # The refmac-refined pdb
        self.name = None
        self.program = None
        self.solution = None
        self.rfact = None
        self.rfree = None
        self.buccRfact = None # Bucc_final_Rfact
        self.buccRfree = None
        self.arpWarpRfact = None # ARP_final_Rfact/Rfree
        self.arpWarpRfree = None
        self.shelxCC = None
        
        self.header = [] # The header format for this table
        
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
    
    def __init__(self, mrbump_dir):
        
        self.mrbump_dir = mrbump_dir
        self.results = []
        
        # List of all the possible column titles and their result object attributes
        self.title2attr = { 
                            'Model_Name' :'name',
                            'MR_Program': 'program',
                            'Solution_Type': 'solution',
                            'final_Rfact' : 'rfact',
                            'final_Rfree' :'rfree',
                            'Bucc_final_Rfact' :'buccRfact',
                            'Bucc_final_Rfree' :'buccRfree',
                            'ARP_final_Rfact' : 'arpWarpRfact',
                            'ARP_final_Rfree' :'arpWarpRfree',
                            'SHELXE_CC' : 'shelxCC' 
                             }
        
        self.logger = logging.getLogger()
        
    def _addFailedResults(self, failed, header):
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
            result.jobDir = os.path.join( self.mrbump_dir, 'search_'+ensemble+'_mrbump' )
            result.header = header
            result.solution = reason
            self._getUnfinishedResult( result )
            self.results.append( result )
            count += 1
        
        self.logger.debug("Added {0} MRBUMP result failures".format(count) )
        return

    def extractResults( self ):
        """
        Find the results from running MRBUMP and sort them
        """

        # Get a list of the ensembles (could get this from the amopt dictionary)
        # For now we just use the submission scripts and assume all have .sub extension
        ensembles = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( self.mrbump_dir, "*.sub") ) ]

        if not len(ensembles):
            self.logger.warn("Could not extract any results from directory: {0}".format( self.mrbump_dir ) )
            return False
        
        # reset any results
        self.results = []
        failed = {} # dict mapping failures to what went wrong - need to process at the end
        header = None
        nfields=None
        for ensemble in ensembles:

            # Check job directory
            jobDir = os.path.join( self.mrbump_dir, 'search_'+ensemble+'_mrbump' )
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
            
            firstLine = True
            
            # Read results table to get the results
            for line in open(resultsTable):
                
                # Create a result object for each line in the output file
                result = MrBumpResult()
                result.jobDir = jobDir
                
                line = line.strip()
                if firstLine:
                    
                    # Processing header
                    firstLine=False
                    header = line.split()
                    nfields = len(header) # count as check
                    badHeader=False
                    for f in header:
                        # Map the data fields to their titles
                        if f not in self.title2attr.keys():
                            result = MrBumpResult()
                            result.jobDir = jobDir
                            self.logger.critical("jobDir {0}: Problem with field {1} in headerline: {2}".format( jobDir, f, line ) )
                            result.header = header
                            result.solution = "problem-header-resultsTable.dat"
                            self._getUnfinishedResult( result )
                            self.results.append( result )
                            badHeader=True
                            break # break out of header check
                    
                    if badHeader:
                        # break out of reading this results table
                        break
                    else:
                        # Got a valid header so continue reading file
                        continue
                    # End header processing
                   
                # Create a result object for each line in the output file
                result = MrBumpResult()
                result.jobDir = jobDir
                    
                result.header = header
                
                fields = line.split()
                if len(fields) != nfields:
                    msg = "jobDir {0}: Problem getting dataline: {1}".format( jobDir, line )
                    self.logger.debug(msg)
                    result.solution = "corrupted-data-resultsTable.dat"
                    self._getUnfinishedResult( result )
                    self.results.append( result )
                    break
                
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
                resultDir = os.path.join( result.jobDir,'data',dirName,'unmod','mr',result.program,'refine' )
                #print resultDir
                result.resultDir = resultDir
                
                # Need to reconstruct final pdb from directory and program name
                pdbName = "refmac_" + result.program + "_loc0_ALL_"  + ensemble + "_UNMOD.pdb"
                result.pdb = os.path.join( resultDir, pdbName )
                
                self.results.append( result )

        if not header or not len(header):
            self.logger.warn("Could not extract any results from directory - no header: {0}".format( self.mrbump_dir ) )
            return False

        # Process the failed results
        if failed:
            self._addFailedResults( failed, header )
                
        if not len(self.results):
            self.logger.warn("Could not extract any results from directory: {0}".format( self.mrbump_dir ) )
            return False

        # Sort the results
        self.sortResults()
        
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
                result.name = "loc0_ALL_" + dname + "_UNMOD"
                result.solution = "ERROR"
            else:
                # Use dirname but remove "loc0_ALL_" from front
                #result.name = os.path.basename( dlist[0] )[9:]
                # Use dirname but add "_UNMOD" to back
                result.name = os.path.basename( dlist[0] )+"_UNMOD"

        result.program = "unknown"
        result.rfact = -1
        result.rfree = -1
        result.shelxCC = -1
        self.buccRfact = -1
        self.buccRfree = -1
        self.arpWarpRfact = -1
        self.arpWarpRfree = -1
        
    def sortResults( self ):
        """
        Sort the results
        """
        
        # Use the first result to determine what attributes are present and how we will sort the results
        reverse=False
        r = self.results[0]
        if r.shelxCC:
            reverse=True
            sortf = lambda x: float( x.shelxCC )
        elif r.buccRfree:
            sortf = lambda x: float( x.buccRfree )
        elif r.arpWarpRfree:
            sortf = lambda x: float( x.arpWarpRfree )
        else:
            sortf = lambda x: float( x.rfree )
        
        # Now sort them
        self.results.sort(key=sortf)
        if reverse:
            self.results.reverse()
    
    def summariseResults( self ):
        """Return a string summarising the results"""
        
        got = self.extractResults()
        if got:
            return self.summaryString()
        else:
            return "\n!!! No results found in directory: {0}\n".format( self.mrbump_dir )
    
    def summaryString( self ):
        """Return a string suitable for printing the sorted results"""
        
        resultsTable = []
        
        #Header
        resultsTable.append( self.results[0].header )
        
        for result in self.results:
            resultLine = []
            for h in result.header:
                resultLine.append( getattr( result, self.title2attr[h] ) )
            resultsTable.append( resultLine )

        # Format the results
        table = printTable.Table()
        summary = table.pprint_table( resultsTable )
        
        r = "\n\nOverall Summary:\n\n"
        r += summary
        r += '\nBest results so far are in :\n\n'
        r +=  self.results[0].resultDir
            
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

    r = ResultsSummary( mrbump_dir )
    print r.summariseResults()
    
