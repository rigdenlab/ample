#!/usr/bin/env python

import glob
import locale
import logging
import os
import re

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
        self.shelxCC = None
        
        self.header = None # The header format for this table

class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """
    
    def __init__(self, mrbump_dir, cluster=False ):
        
        self.mrbump_dir = mrbump_dir
        self.cluster = cluster # whether jobs were run on a cluster
        self.results = []
        
        self.logger = logging.getLogger()

    def extractResults( self ):
        """
        Find the results from running MRBUMP and sort them
        """

        # reset any results
        self.results = []
        
        # how we recognise a job directory
        dir_re = re.compile("^search_.*_mrbump$")
        
        jobDirs = []
        if self.cluster:
            for pd in os.listdir( self.mrbump_dir ):
                pd = os.path.join( self.mrbump_dir, pd )
                if os.path.isdir( pd ):
                    for d in os.listdir( pd ):
                        dpath = os.path.join( self.mrbump_dir, pd, d )
                        if dir_re.match( d ) and os.path.isdir( dpath ):
                            jobDirs.append( dpath )
        else:
            jobDirs = glob.glob( os.path.join( self.mrbump_dir, "search_*_mrbump" ) )
        
        if not len(jobDirs):
            self.logger.warn("Could not extract any results from directory: {0}".format( self.mrbump_dir ) )
            return False
        
        header = None
        for jobDir in jobDirs:
            
            self.logger.debug(" -- checking directory for results: {0}".format( jobDir ) )
            resultsTable = os.path.join( jobDir,"results", "resultsTable.dat" )
            if not os.path.exists(resultsTable):
                self.logger.debug(" -- Could not find file: {0}".format( resultsTable ) )
                continue
            
            firstLine = True
            # Read results table to get the results
            for line in open(resultsTable):
                
                line = line.strip()
                
                if firstLine:
                    # probably overkill...
                    if line != "Model_Name   MR_Program   Solution_Type   final_Rfact   final_Rfree   SHELXE_CC":
                        raise RuntimeError,"Problem getting headerline: {0}".format(line)
                    header = line
                    firstLine=False
                    continue
                
                result = MrBumpResult()
                result.jobDir = jobDir
                
                fields = line.split()
                
                result.name = fields[0]
                result.program = fields[1].lower()
                result.solution = fields[2]
                result.rfact = fields[3]
                result.rfree = fields[4]
                result.shelxCC = fields[5]
                result.header = header
                
                if result.program not in ['phaser','molrep']:
                    raise RuntimeError,"getResult, unrecognised program in line: {0}".format(line)
                
                # Rebuild the path that generated the result
                # Strip  _UNMOD from (e.g.): loc0_ALL_All_atom_trunc_0.34524_rad_1_UNMOD 
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
 
    def resultsTable(self, table):
        """Returns a string of table data, padded for alignment
        @param table: The table to print. A list of lists.
        Each row must have the same number of columns. 
        
        """
        
        col_paddings = []
        
        out = ""
    
        for i in range(len(table[0])):
            col_paddings.append(self._get_max_width(table, i))
    
        for row in table:
            # left col
            out += row[0].ljust(col_paddings[0] + 1)
            # rest of the cols
            for i in range(1, len(row)):
                col = self._format_num(row[i]).rjust(col_paddings[i] + 2)
                out +=  col
            out += "\n"
                
        return out
    
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
        resultsTable.append( self.results[0].header.split() )
        
        for result in self.results:
            rl = [ result.name,
                    result.program,
                    result.solution,
                    result.rfact,
                    result.rfree,
                    result.shelxCC,
                  ]
            resultsTable.append( rl )
    
        # Format the results
        summary = self.resultsTable( resultsTable )
        
        r = "\n\nOverall Summary:\n\n"
        r += summary
        r += '\nBest results so far are in :\n\n'
        r +=  self.results[0].resultDir
            
        return r

    def _format_num(self, num):
        """Format a number according to given places.
        Adds commas, etc. Will truncate floats into ints!"""
    
        try:
            if "." in num:
                inum = float(num)
                return locale.format("%.2f", (0, inum), True)
            else:
                inum = int(num)
                return locale.format("%.*f", (0, inum), True)
    
        except (ValueError, TypeError):
            return str(num)
    
    def _get_max_width(self, table, index):
        """Get the maximum width of the given column index"""
        return max([len(self._format_num(row[index])) for row in table])
    
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
    if len(sys.argv) == 3:
        cluster = True
        mrbump_dir = os.path.join( os.getcwd(), sys.argv[2] )
    elif len(sys.argv) == 2:
        mrbump_dir = os.path.join( os.getcwd(), sys.argv[1] )
        cluster=False
    else:
        mrbump_dir = "/Users/jmht/Documents/AMPLE/res.test/cluster_run1"
        mrbump_dir = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/MRBUMP/cluster_1"
        mrbump_dir = "/gpfs/home/HCEA041/djr01/jxt15-djr01/TM/3OUF/ROSETTA_MR_1/MRBUMP/cluster_run1"
    
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    r = ResultsSummary( mrbump_dir, cluster=cluster )
    print r.summariseResults()
    
