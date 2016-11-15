
class ClusterData(object):
    """Class to hold all data related to a cluster"""
    def __init__(self):
        self.cluster_method = None
        self.index = None # index of the cluster in the list of clusters (counting from 1)
        self.models = []
        self.num_clusters = None # How many clusters there are in list of clusters
        self.r_cen = []  # ordered list of the distance from the cluster centroid for each pdb (for Spicker)
        self.score_type = None
        
        # This may be set if it's not the first model
        self._centroid = None
    
    @property
    def centroid(self):
        return self.models[0] if self._centroid is None else self._centroid
    
    @centroid.setter
    def centroid(self, model):
        self._centroid = model
    
    @property
    def size(self):
        return self.__len__()
    
    def __len__(self):
        """Return the number of models in this cluster."""
        return 0 if self.models is None else len(self.models)

    def __str__(self):
        """Return a string representation of this object."""
        _str = super(Cluster, self).__str__() + "\n"
        # Iterate through all attributes in order
        for k in sorted(self.__dict__.keys()):
            _str += "{0} : {1}\n".format(k, self.__dict__[k])
        return _str
