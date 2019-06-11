import os

class RosettaScoreData(object):
    
    def __init__(self):
        self.score = None
        self.rms = None
        self.maxsub = None
        self.description = None
        self.model = None


class RosettaScoreParser(object):
    
    def __init__(self, directory ):
        
        self.directory = directory
        
        self.avgScore = None
        self.topScore = None
        self.avgRms = None
        self.topRms = None
        self.avgMaxsub = None
        self.topMaxsub = None
        self.data = []
        
        score_file = os.path.join( directory, "score.fsc")
        if not os.path.isfile(score_file):
            raise RuntimeError("Cannot find ROSETTA score file: {0}".format(score_file))
        self.parseFile(score_file)
        
    def parse_file(self, score_file):
        idxScore=None
        idxRms=None
        idxMaxsub=None
        idxDesc=None
        for i, line in enumerate(open(score_file, 'r')):
            line = line.strip()
            # Read header
            if i == 0:
                for j,f in enumerate(line.split()):
                    if f=="score":
                        idxScore=j
                    elif f=="rms":
                        idxRms=j
                    elif f=="maxsub":
                        idxMaxsub=j
                    elif f=="description":
                        idxDesc=j
                
                if idxScore is None or idxRms is None or idxMaxsub is None or idxDesc is None:
                    raise RuntimeError("Missing header field from score file: {0}".format(score_file))
                continue
                # End read header
    
            if not line: # ignore blank lines - not sure why they are there...
                continue
            
            d = RosettaScoreData()
            
            fields = line.split()
            d.score = float(fields[idxScore])
            d.rms = float(fields[idxRms])
            d.maxsub = float(fields[idxMaxsub])
            d.description = fields[idxDesc]
            #pdb = fields[31]
            
            d.model = os.path.join(self.directory, d.description+".pdb")
            
            self.data.append(d)
        
        avg = 0
        self.topScore = self.data[0].score
        for d in self.data:
            avg += d.score
            if d.score < self.topScore:
                self.topScore = d.score
        self.avgScore  = avg / len(self.data)
        
        avg = 0
        self.topRms = self.data[0].rms
        for d in self.data:
            avg += d.rms
            if d.rms < self.topRms:
                self.topRms = d.rms
        self.avgRms  = avg / len(self.data)
        
        avg = 0
        self.topMaxsub = self.data[0].maxsub
        for d in self.data:
            avg += d.maxsub
            if d.maxsub > self.topMaxsub:
                self.topMaxsub = d.maxsub
        self.avgMaxsub  = avg / len(self.data)
        
        return
        
    def maxsub_sorted(self, reverse=True):
        return sorted(self.data, key=lambda data: data.maxsub, reverse=reverse)
     
    def rms_sorted(self, reverse=True):
        return sorted(self.data, key=lambda data: data.rms, reverse=reverse)
    
    def rms(self, name):
        for d in self.data:
            if d.description == name:
                return d.rms
            
    def maxsub(self, name):
        for d in self.data:
            if d.description == name:
                return d.maxsub
    
    def __str__(self):
        s = "Results for: {0}\n".format(self.name)
        s += "Top score : {0}\n".format(self.topScore)
        s += "Avg score : {0}\n".format(self.avgScore)
        s += "Top rms   : {0}\n".format(self.topRms)
        s += "Avg rms   : {0}\n".format(self.avgRms)
        s += "Top maxsub: {0}\n".format(self.topMaxsub)
        s += "Avg maxsub: {0}\n".format(self.avgMaxsub)
        return s
