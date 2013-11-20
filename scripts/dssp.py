import unittest

class DsspParser(object):
    """
    Class 
    """

    def __init__(self,pfile):

        self.dsspfile = pfile
        
        self.chainIds = []
        self.resNames = []
        self.resSeqs = []
        self.assignment = []
        
        self.percentH = []
        self.percentC = []
        self.percentE = []

        self.parse()

        return
    
    def parse(self):
        """Info from: http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html
        """

        # print os.path.join(os.getcwd(), logfile)

        self.chainIds = []
        self.resNames = []
        self.resSeqs = []
        self.assignment = []
        
        capture=False
        currentChain = None
        for line in open(self.dsspfile, 'r'):
            
            if "#  RESIDUE" in line:
                capture=True
                continue
                
            if capture:
                
                # Ignore chain break characters - we use the chainId
                if "!" in line:
                    continue
                
                #print "\"{0}\"".format(line)
                #idx = int( line[0:5].strip() )
                resSeq = int( line[5:10].strip() )
                chainId = line[10:12].strip()
                resName = line[12:14].strip()
                #print "\"{0}\"".format(line[14:17])
                assign = line[16]
                 
                if currentChain != chainId:
                    currentChain = chainId
                    self.chainIds.append( chainId )
                    self.resNames.append( [] )
                    self.resSeqs.append( [] )
                    self.assignment.append( [] )
                                        
                self.resNames[-1].append( resName )
                self.resSeqs[-1].append( resSeq )
                self.assignment[-1].append( assign )
                
        if not len( self.resNames[0] ) or not len( self.assignment[0] ):
            raise RuntimeError,"Got no assignment!"
         
        for chain in range( len( self.chainIds ) ):
            nH = 0
            nC = 0
            nE = 0
            for p in self.assignment[chain]:
                if p == "H":
                    nH += 1
                elif p == "E":
                    nE += 1
                # Just assume everything else is a coil
                else:
                    nC += 1
            
            self.percentC.append( float(nC) / len(self.assignment[ chain ] ) * 100 )
            self.percentH.append( float(nH) / len(self.assignment[ chain ] ) * 100 )
            self.percentE.append( float(nE) / len(self.assignment[ chain ] ) * 100 )
        
        return
    
    def asDict(self):
        d = {}
        d['chainIds'] = self.chainIds
        d['assignment'] = self.assignment
        d['resNames'] = self.resNames
        d['resSeqs'] = self.resSeqs
        d['percentC'] = self.percentC
        d['percentE'] = self.percentE
        d['percentH'] = self.percentH
        
        return d
    
    def getAssignment(self, resSeq, resName, chainId ):
        ci = self.chainIds.index( chainId )
        ri = self.resSeqs[ ci ].index( resSeq )
        # Just a check to make sure things are working - ignore X as it'll be a non-standard residue e.g. N-FORMYLMETHIONINE
        dsspResName = self.resNames[ ci ][ ri ]
        # in dssp cysteine bridges are signified by lower-case letters
        if dsspResName != resName and not dsspResName.islower() and dsspResName != 'X' :
            raise RuntimeError,"Missmatching residues: {0}: {1}".format( self.resNames[ ci ][ ri ], resName )
        return self.assignment[ ci ][ ri ]


class Test(unittest.TestCase):

    def testParse1(self):
        """parse 2bhw"""
        
        dssp_file = "/media/data/shared/TM/2BHW/2bhw.dssp"
        dsspP = DsspParser( dssp_file )
        a = [' ', ' ', 'T', 'T', ' ', 'T', 'T', 'S', 'S', 'T', 'T', ' ', ' ', ' ', 'T', 'T', 'G', 'G', 'G', ' ', ' ', 'S', ' ', ' ', 'T', 'T', ' ', ' ', 'S', ' ', 'S', 'T', 'T', ' ', ' ', 'S', ' ', ' ', 'T', 'T', ' ', 'T', 'T', ' ', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'T', 'T', 'T', ' ', ' ', ' ', 'S', ' ', ' ', 'S', 'G', 'G', 'G', 'S', 'G', 'G', 'G', 'G', 'G', 'S', 'T', 'T', ' ', 'E', 'E', 'G', 'G', 'G', ' ', 'T', 'T', 'S', 'E', 'E', 'E', ' ', ' ', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'E', 'E', 'T', 'T', 'E', 'E', ' ', 'S', 'S', 'S', 'S', ' ', ' ', ' ', 'T', 'T', 'S', ' ', 'T', 'T', ' ', 'T', 'T', ' ', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', ' ', 'S', ' ', 'S', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', ' ', 'T', 'T', 'T', 'S', 'S', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', ' ', 'S', ' ', ' ']
        self.assertEqual( dsspP.assignment[0], a)
        
        return
        
    def testParse2(self):
        dssp_file = "/media/data/shared/TM/3OUF/3ouf.dssp"
        dsspP = DsspParser( dssp_file )
        
        return
    
if __name__ == "__main__":
    
    unittest.main()


