#!/usr/bin/env python

import glob
import locale
import os
import re

import printTable


def get_rfree(pdb):
  pdb = open(pdb)

  freer_pattern = re.compile('REMARK\s*\d*\s*FREE R VALUE\s*\:\s*(\d*\.\d*)')
  FreeR = 'nan' 

  for line in pdb:

    if re.search(freer_pattern, line):
    #  print line
      split=  re.split(freer_pattern, line)
   #   print split
      FreeR = split[1]

  return FreeR
###########################

def get_shelx_score(RESULT):

    Best_result = 0
    result = open(RESULT)
    fail_pattern = ('\*\* Unable to trace map - giving up \*\*')
    pattern = re.compile('CC for partial structure against native data =\s*(\d*\.\d*)\s*%')
    for line in result: 
       if re.search(pattern, line):
        split= re.split(pattern, line)
      #  print split
        if float(split[1]) > Best_result:
          Best_result = float(split[1])
       if re.search(fail_pattern, line):
        Best_result = line.rstrip('\n')

    return Best_result
  
##########################
def make_log(mr_bump_path, final_log):
  final_log = open(final_log, "w")

  list_dir = os.listdir(mr_bump_path)
  for a_dir in list_dir:
    if os.path.isdir(mr_bump_path + '/'+a_dir):

     name=re.sub('search_', '', a_dir)
     name=re.sub('_mrbump', '', name)
     
     phaser_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'
     molrep_pdb = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'
     phaser_log = mr_bump_path + '/'+a_dir+'/data/loc0_ALL_'+name+'/unmod/mr/phaser/phaser_loc0_ALL_'+name+'_UNMOD.1.pdb'

     shelx_phaser = mr_bump_path + '/'+a_dir+'/phaser_shelx/RESULT'
     shelx_molrep = mr_bump_path + '/'+a_dir+'/molrep_shelx/RESULT'

     if os.path.exists(phaser_pdb):
      if os.path.exists(shelx_phaser):
        phaser_FreeR =   get_rfree(phaser_pdb)  
        phaser_shelx =   get_shelx_score(shelx_phaser)
        final_log.write('Ensembe ' + name+'  phaser FreeR: '+phaser_FreeR +  '  shelx score: '+str(phaser_shelx) + '\n')
        final_log.flush()
  
     if os.path.exists(molrep_pdb):
      if os.path.exists(shelx_molrep):
        molrep_FreeR =   get_rfree(molrep_pdb) 
        molrep_shelx =   get_shelx_score(shelx_molrep)  
        final_log.write('Ensembe ' + name+'  molrep FreeR: '+molrep_FreeR +  '  shelx score: '+str(molrep_shelx) + '\n')
        final_log.flush()

############
#mr_bump_path = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP'
#final_log = '/home/jaclyn/Baker/test_cases/1MB1/MRBUMP/FINAL'
#make_log(mr_bump_path,final_log)


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
    
    def __init__(self):
        
        self.mrbump_dir = None
        self.cluster = None # whether jobs were run on a cluster
        self.results = []

    def summariseResults(self, mrbump_dir, cluster=False ):
        """Return a string summarising the results"""
        
        self.mrbump_dir = mrbump_dir
        self.cluster = cluster
        
        self.extractResults()
        self.sortResults()
        
        return self.summaryString()

    def extractResults( self ):
        """
        Find the results from running MRBUMP
        """
        
        # how we recognise a job directory
        dir_re = re.compile("^search_.*_mrbump$")
        
        jobDirs = []
        if self.cluster:
            for pd in os.listdir( self.mrbump_dir ):
                if os.path.isdir( pd ):
                    for d in os.listdir( os.path.join(self.mrbump_dir, pd) ):
                        if dir_re.match( d ) and os.path.isdir( d ):
                            jobDirs.append( os.path.join( self.mrbump_dir, pd, d )  )
        else:
            jobDirs = glob.glob( os.path.join( self.mrbump_dir, "search_*_mrbump" ) )
        
        if not len(jobDirs):
            raise RuntimeError,"Could not extract any results"
        
        header = None
        for jobDir in jobDirs:
            
            resultsTable = os.path.join( jobDir,"results", "resultsTable.dat" )
            if not os.path.exists(resultsTable):
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
            raise RuntimeError,"Could not extract any results"
    
    def sortResults( self ):
        """
        Sort the results and return an ordered list
        """
        
        use_shelx=True # For time being assume we always use shelx
        SHELXSUCCESS = 25.0
        
        # sort results
        # Use header of first result to work out index
        if use_shelx:
            sortf = lambda x: float( x.shelxCC )
        else:
            sortf = lambda x: float( x.rfree )
        
        self.results.sort(key=sortf)
        self.results.reverse()
    
    def summaryString( self ):
        """
        Return a string suitable for printing the sorted results"""
        
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
        r += '\nBest results so far are in :\n'
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


if __name__ == "__main__":
    mrbump_dir = "/Users/jmht/Documents/AMPLE/res.test/cluster_run1"
    mrbump_dir = "/Users/jmht/Documents/AMPLE/res.test/MRBUMP"
    
    r = ResultsSummary()
    print r.summariseResults( mrbump_dir, cluster=False )
    