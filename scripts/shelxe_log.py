'''
Created on 28 Nov 2013

@author: jmht
'''

import os
import re

class ShelxeLogParser(object):
    """
    Class to mine information from a shelxe log
    """
    
    def __init__(self,logfile):
        
        self.logfile = logfile
        
        # All data refers to the best cycle
        self.CC = None
        self.avgChainLength = None
        self.maxChainLength = None
        self.numChains = None
        self.cycle = None
        
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
                (cc, avgChainLength, maxChainLength, numChains ) = self._parseCycle(fh)
                cycleData.append( (cc, avgChainLength, maxChainLength, numChains) )
            
            
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
                self.maxChainLength = cycleData[ cycle-1 ][2]
                self.numChains = cycleData[ cycle-1 ][3]
                self.cycle = cycle

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
        maxChainLength = max(lengths)  
        numChains = len(lengths)    
        
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
        
        return ( cc, avgChainLength, maxChainLength, numChains )

#END ShelxeLogParser
