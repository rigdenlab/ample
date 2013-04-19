'''
Created on Apr 18, 2013

@author: jmht
'''

import logging
import os

class EnsembleData(object):
    """Class to hold data about an ensemble"""
    
    def __init(self):
        
        self.name = None 
        self.numModels = None
        self.numResidues = None
        self.numAtoms = None
        self.residues = []
        
        self.sideChainTreatment = None
        self.radiusThreshold = None
        self.truncationThreshold = None
        
        self.pdb = None # path to the ensemble file


class Ensembler(object):
    """Class to generate ensembles from cluster of models"""
    
    def __init__(self):
        """Initialise"""
        
        self.workdir = None
        self.percent= None
        # List of the thresholds we truncate at
        self.truncation_thresholds = []
        # List of directories where the truncated files are (must match truncation_thresholds )
        self.truncation_dirs = []
    
    def genereate_ensembles(self, cluster_models=None, workdir=None, percent=None ):
        """Generate the ensembles for this list of cluster models"""
        
        if not os.path.exists( workdir ):
            os.mkdir( workdir )
        self.workdir = workdir
        os.chdir( self.workdir )
        
        self.truncate_models( percent=percent )
        
        # Make sure we got something to work with
        assert len( self.truncation_thresholds ) and len( self.truncation_thresholds ) == len( self.truncation_dirs ), \
        "Must have some thresholds and directories for each truncation level!"
        
        for i, threshold in enumerate( self.truncation_thresholds ):
            self.make_truncated_ensemble( threshold, self.truncation_dirs[i] )
            
        self.edit_sidechains()
        
        
    def truncate_models(self, percent=None):
        """Truncate the models according to the percent criteria"""
        pass
        
    def make_truncated_ensemble(self, truncation_threshold, truncation_dir ):
        pass
    
    def edit_sidechains(self):
        """Take the ensembles and give them the 3 sidechain treatments"""
        pass
    

if __name__ == "__main__":
    
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    cf="/Users/jmht/Documents/AMPLE/ample-dev1/examples/toxd-example/ROSETTA_MR_0/spicker_run/spicker_cluster_1.list"
    cluster_models = [ m for m in open( cf, "r" ) ]
    workdir="/Users/jmht/Documents/AMPLE/ample-dev1/examples/toxd-example/jtest"
    percent=10
    
    ensembler = Ensembler()
    ensembler.generate_ensembles( cluster_models=cluster_models, workdir=workdir, percent=percent )
    ensembles = ensembler.get_file_list()
    print ensembles



