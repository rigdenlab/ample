#!/usr/bin/env python

import glob
import locale
import logging
import os
import re

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

class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """
    
    def __init__(self, mrbump_dir):
        
        self.mrbump_dir = mrbump_dir
        self.results = []
        
        # List of all the possible column titles
        self.columnTitles = [ 'Model_Name',
                             'MR_Program',
                             'Solution_Type',
                             'final_Rfact',
                             'final_Rfree',
                             'Bucc_final_Rfact',
                             'Bucc_final_Rfree',
                             'ARP_final_Rfact',
                             'ARP_final_Rfree',
                             'SHELXE_CC' ]
        
        # list of attributes of the result object that match the columnTitles (same order)
        self.resultAttr = [ 'name',
                           'program',
                           'solution',
                           'rfact',
                           'rfree',
                           'buccRfact',
                           'buccRfree',
                           'arpWarpRfact',
                           'arpWarpRfree',
                           'shelxCC' ]
        
        self.logger = logging.getLogger()

    def extractResults( self ):
        """
        Find the results from running MRBUMP and sort them
        """

        # how we recognise a job directory
        dir_re = re.compile("^search_.*_mrbump$")
        jobDirs = []
        #jobDirs = glob.glob( os.path.join( self.mrbump_dir, "search_*_mrbump" ) )
        for adir in os.listdir( self.mrbump_dir ):
            # REM only returns relative path
            dpath = os.path.join( self.mrbump_dir, adir )
            if dir_re.match( adir ) and os.path.isdir( dpath ):
                jobDirs.append( dpath )
        
        if not len(jobDirs):
            self.logger.warn("Could not extract any results from directory: {0}".format( self.mrbump_dir ) )
            return False
        
        # reset any results
        self.results = []
        header = None
        nfields=None
        for jobDir in jobDirs:
            
            self.logger.debug(" -- checking directory for results: {0}".format( jobDir ) )
                        
            # Check if finished
            if not os.path.exists( os.path.join( jobDir, "results", "finished.txt" ) ):
                self.logger.debug(" Found unfinished job: {0}".format( jobDir ) )
                if header:
                    result = MrBumpResult()
                    result.jobDir = jobDir
                    result.solution = "unfinished"
                    result.header = header
                    self.getUnfinishedResult( result )
                    self.results.append( result )
                else:
                    self.logger.critical(" Cannot provide result for: {0}".format( jobDir ) )
                continue
            
            resultsTable = os.path.join( jobDir,"results", "resultsTable.dat" )
            if not os.path.exists(resultsTable):
                self.logger.debug(" -- Could not find file: {0}".format( resultsTable ) )
                if header:
                    result = MrBumpResult()
                    result.jobDir = jobDir
                    result.header = header
                    result.solution = "no-resultsTable.dat"
                    self.getUnfinishedResult( result )
                    self.results.append( result )
                else:
                    self.logger.critical(" Cannot provide result for: {0}".format( jobDir ) )
                continue
            
            firstLine = True
            # This maps the index of data field to the index of the columnTitle and resultAttr 
            fieldIndex = [ None ] * len( self.columnTitles )
            # Read results table to get the results
            for line in open(resultsTable):
                line = line.strip()
                if firstLine:
                    firstLine=False
                        
                    # Processing header
                    header = line.split()
                    nfields = len(header) # count as check
                    herror=None
                    for i, f in enumerate( header ):
                        # Map the data fields to their titles
                        try:
                            findex = self.columnTitles.index( f )
                            fieldIndex[ findex ] = i
                        except ValueError:
                            self.logger.critical("jobDir {0}: Problem getting headerline: {1}".format( jobDir, line ) )
                            result.header = header
                            result.solution = "corrupted-header-resultsTable.dat"
                            self.getUnfinishedResult( result )
                            self.results.append( result )
                            break
                    continue
                    # End header processing
                
                # Create result object
                result = MrBumpResult()
                result.header = header
                result.jobDir = jobDir
                
                fields = line.split()
                if len(fields) != nfields:
                    msg = "jobDir {0}: Problem getting dataline: {1}".format( jobDir, line )
                    result.solution = "corrupted-data-resultsTable.dat"
                    self.getUnfinishedResult( result )
                    self.results.append( result )
                    break
                
                # Now set all attributes of the result object
                for i, f in enumerate( fields ):
                    # Strip loc0_ALL_ from front and strip  _UNMOD from end from (e.g.): loc0_ALL_All_atom_trunc_0.34524_rad_1_UNMOD
                    #result.name = fields[0][9:-6]
                    # Don't do the above yet till we've finsihed the next set of runs
                    findex = fieldIndex.index( i )
                    attr = self.resultAttr[findex]
                    if attr == 'program':
                        setattr( result, attr, f.lower() )
                    else:
                        setattr( result, attr, f )

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
                self.results.append( result )
                
        if not len(self.results):
            self.logger.warn("Could not extract any results from directory: {0}".format( self.mrbump_dir ) )
            return False
        
        # Sort the results
        self.sortResults()
        
        return True
    
    def getUnfinishedResult(self, result ):
        """Return a result for an unfinished job"""
        
        
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
        self.buccRfact = -1
        self.buccRfree = -1 
        
    def sortResults( self ):
        """
        Sort the results
        """
        
        use_shelx=True # For time being assume we always use shelx
        if use_shelx:
            sortf = lambda x: float( x.shelxCC )
        else:
            sortf = lambda x: float( x.rfree )
        
        self.results.sort(key=sortf)
        self.results.reverse()
    
    def summariseResults( self ):
        """Return a string summarising the results"""
        
        got = self.extractResults()
        if got:
            return self.summaryString()
        else:
            return "\n!!! No results found in directory: {0}\n".format( mrbump_dir )
    
    def summaryString( self ):
        """Return a string suitable for printing the sorted results"""
        
        resultsTable = []
        
        #Header
        resultsTable.append( self.results[0].header )
        
        for result in self.results:
            resultLine = []
            for h in result.header:
                i = self.columnTitles.index( h )
                resultLine.append( getattr( result, self.resultAttr[i] ) )
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
    
