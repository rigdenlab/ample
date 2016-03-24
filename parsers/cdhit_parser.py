
class CDhitLogParser(object):
    """ Class to mine information from a cdhit clstr file """

    def __init__(self):
        self.clusters = {}
        self.nrClusters = 0
        self.nrEffSeqs = 0

    def parse(self, clstrfile):
        self.clusters = self.getClusters(clstrfile)

        self.nrClusters = len(self.clusters.keys())
        self.nrEffSeqs = int(sum(i['neffWeight'] for i in self.clusters.values()))
        return

    def getClusters(self, clstrfile):
        clusters = {}
        
        with open(clstrfile, 'r') as fh:
            line = fh.readline()
            while True:
                if not line: break
                
                if line.startswith(">"):
                    cluster_idx = int(line.strip().split()[-1])
                    clusters[cluster_idx] = {'count': 0, 'averageId': 100.00, 'averageLength': 0, 'centroid': None, 'neffWeight': 1.00} 
                else:
                    clusters[cluster_idx]['count'] += 1

                    line = line.strip().split()
                    seqId = line[2][1:]
                    identity = line[-1][:-1]
                    length   = line[1][:-3]
                    
                    if identity and "*" not in identity:
                        clusters[cluster_idx]['averageId'] += float(identity)
                    if identity and "*" in identity:
                        clusters[cluster_idx]['centroid'] += seqId

                    clusters[cluster_idx]['averageLength'] += float(length)

                line = fh.readline()
       
        return self.averageClusters(clusters)
        
    def averageClusters(self, clusters):
        for cluster_data in clusters.values():
            cluster_data['averageId'] = cluster_data['averageId'] / cluster_data['count']
            cluster_data['averageLength'] = abs(float(cluster_data['averageLength']) / cluster_data['count'])
            cluster_data['neffWeight'] = cluster_data['neffWeight'] / cluster_data['count']
        return clusters
