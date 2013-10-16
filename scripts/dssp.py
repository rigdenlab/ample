

class DsspParser(object):
    """
    Class 
    """

    def __init__(self,pfile):

        self.dsspfile = pfile
        
        self.resNames = []
        self.assignment = []
        
        self.percentH = None
        self.percentC = None
        self.percentE = None

        self.parse()

        return
    
    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)

        self.resNames = []
        self.resSeqs = []
        self.assignment = []
        
        capture=False
        for line in open(self.dsspfile, 'r'):
            
            if "#  RESIDUE" in line:
                capture=True
                continue
                
            if capture:
                
                # Only capture first chain
                if "!" in line:
                    capture=False
                    break
                
                #print "\"{0}\"".format(line)
                #idx = int( line[0:5].strip() )
                resSeq = int( line[5:10].strip() )
                #chainId = line[10:12].strip()
                resName = line[12:14].strip()
                #print "\"{0}\"".format(line[14:17])
                assign = line[16]
                 
                self.resNames.append( resName )
                self.resSeqs.append( resSeq )
                self.assignment.append( assign )
                
        if not len(self.resNames) or not len( self.assignment):
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
        d['resNames'] = self.resNames
        d['percentC'] = self.percentC
        d['percentE'] = self.percentE
        d['percentH'] = self.percentH
        
        return d
    
    
if __name__ == "__main__":
    dssp_file = "/media/data/shared/TM/2BHW/2bhw.dssp"
    dsspP = DsspParser( dssp_file )
    print dsspP.asDict()