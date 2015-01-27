'''
Created on Apr 18, 2013

This structure of the ensembling modules is dictated by the need to be able to pickle
and unpickle the results dictionary. As such all objects need to have a qualified path
(e.g. ample_ensemble.Ensemble ) - otherwise, the module is taken as main, so that when
the results are unpickled, it will look for Ensemble as an attribute of the main ample
script.

@author: jmht
'''

import cPickle
import glob
import logging
import os
import sys
import unittest

# our imports
import ample_ensemble
import ample_util
import printTable

def create_ensembles( amoptd ):
    """Create the ensembles using the values in the amoptd dictionary"""

    ensembler = ample_ensemble.Ensembler()
    
    ensembler.theseus_exe = amoptd['theseus_exe'] 
    ensembler.maxcluster_exe = amoptd['maxcluster_exe'] 
    ensembler.subcluster_exe = amoptd['maxcluster_exe']
    ensembler.max_ensemble_models = amoptd['max_ensemble_models']
    if amoptd['cluster_method'] == 'spicker':
        cluster_exe = amoptd['spicker_exe']
    
    work_dir=os.path.join(amoptd['work_dir'],'ensembling')
    os.mkdir(work_dir)
    os.chdir(work_dir)
        
    models=glob.glob(os.path.join(amoptd['models_dir'],"*.pdb"))
    ensembles=ensembler.generate_ensembles(models,
                                           cluster_method=amoptd['cluster_method'] ,
                                           cluster_exe=cluster_exe,
                                           num_clusters=amoptd['num_clusters'] ,
                                           percent_truncation=amoptd['percent'],
                                           truncation_method=amoptd['truncation_method'],
                                           work_dir=work_dir)
    
    amoptd['ensembles'] = ensembles
    amoptd['ensembles_data'] = ensembler.ensembles_data

    return

def ensemble_summary(ensembles_data):
    """Print a summary of the ensembling process"""

    clusters={}
    # Loop through all ensemble data objects and build up a data tree
    cluster_method=None
    truncation_method=None
    percent_truncation=None
    for e in ensembles_data:
        if not cluster_method:
            cluster_method = e.cluster_method
            percent_truncation = e.percent_truncation
            truncation_method = e.truncation_method
            num_clusters = e.num_clusters
        if e.cluster_num not in clusters:
            clusters[e.cluster_num] = {}
            clusters[e.cluster_num]['cluster_centroid'] = e.cluster_centroid
            clusters[e.cluster_num]['cluster_num_models'] = e.cluster_num_models
            clusters[e.cluster_num]['tlevels'] = {}
            
        if e.truncation_level not in clusters[e.cluster_num]['tlevels']:
            clusters[e.cluster_num]['tlevels'][e.truncation_level] = {}
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['truncation_variance']=e.truncation_variance
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['num_residues']=e.num_residues
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds']={}
        
        if e.subcluster_radius_threshold not in clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds']:
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]={}
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['num_models']=e.num_models
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['sct']={}

        if e.side_chain_treatment not in clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['sct']:
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['sct'][e.side_chain_treatment]={}
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['sct'][e.side_chain_treatment]['name']=e.name
            clusters[e.cluster_num]['tlevels'][e.truncation_level]['radius_thresholds'][e.subcluster_radius_threshold]['sct'][e.side_chain_treatment]['num_atoms']=e.num_atoms
    
    tableFormat = printTable.Table()
    rstr = "\n"
    rstr += "Ensemble Results\n"
    rstr += "----------------\n\n"
    rstr += "Cluster method: {0}\n".format(cluster_method)
    rstr += "Truncation method: {0}\n".format(truncation_method)
    rstr += "Percent truncation: {0}\n".format(percent_truncation)
    rstr += "Number of clusters: {0}\n".format(num_clusters)
    
    for cn in sorted(clusters.keys()):
        rstr += "\n"
        rstr += "Cluster {0}\n".format(cn)
        rstr += "Number of models: {0}\n".format(clusters[cn]['cluster_num_models'])
        rstr += "Cluster centroid: {0}\n".format(clusters[cn]['cluster_centroid'])
        rstr += "\n"
        
        tdata = [("Nane","Truncation Level", "Variance Threshold (A^2)", "No. Residues", "Radius Threshold (A)", "No. Decoys", "Number of Atoms", "Sidechain Treatment" ) ]
        for tl in sorted(clusters[cn]['tlevels']):
            tvar=clusters[cn]['tlevels'][tl]['truncation_variance']
            nresidues=clusters[cn]['tlevels'][tl]['num_residues']
            for i, rt in enumerate(sorted(clusters[cn]['tlevels'][tl]['radius_thresholds'])):
                nmodels=clusters[cn]['tlevels'][tl]['radius_thresholds'][rt]['num_models']
                for j, sct in enumerate(sorted(clusters[cn]['tlevels'][tl]['radius_thresholds'][rt]['sct'])):
                    name=clusters[cn]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['name']
                    num_atoms=clusters[cn]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['num_atoms']
                    if i==0:
                        if j==0:
                            tdata.append((name, tl, tvar, nresidues, rt, nmodels, num_atoms, sct))
                        else:    
                            tdata.append((name, "", "", "", "", "", num_atoms, sct))
                    else:
                        tdata.append((name, "", "", "", rt, nmodels, num_atoms, sct))   
                    
        rstr += tableFormat.pprint_table( tdata )        
    
    rstr += "\nGenerated {0} ensembles\n\n".format(len(ensembles_data))
    
    return rstr

class Test(unittest.TestCase):

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ample_dir = os.sep.join( paths[ : -1 ] )
        self.theseus_exe=ample_util.find_exe('theseus')
        self.spicker_exe=ample_util.find_exe('spicker')
        self.maxcluster_exe=ample_util.find_exe('maxcluster')

        return
    
    def testSummary(self):

        ensembler=ample_ensemble.Ensembler()

        work_dir=os.path.join(os.getcwd(),"summary")
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        models=glob.glob(mdir+os.sep+"*.pdb")

        num_clusters=1
        cluster_method='spicker'
        percent_truncation=20
        truncation_method="percent"
        ensembler.generate_ensembles(models,
                                     cluster_method=cluster_method,
                                     cluster_exe=self.spicker_exe,
                                     num_clusters=num_clusters,
                                     percent_truncation=percent_truncation,
                                     truncation_method=truncation_method,
                                     work_dir=work_dir)
        
        print ensemble_summary(ensembler.ensembles_data)
        
        return


if __name__ == "__main__":

    # This runs the ensembling starting from a pickled file containing an amopt dictionary.
    # - used when submitting the modelling jobs to a cluster

    if len(sys.argv) != 2 or not os.path.isfile( sys.argv[1] ):
        print "ensemble script requires the path to a pickled amopt dictionary!"
        sys.exit(1)

    # Get the amopt dictionary
    with open(sys.argv[1], "r") as f:
        amoptd = cPickle.load(f)

    #if os.path.abspath(fpath) != os.path.abspath(amoptd['results_path']):
    #    print "results_path must match the path to the pickle file"
    #    sys.exit(1)

    # Set up logging - could append to an existing log?
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fl = logging.FileHandler("ensemble.log")
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)

    # Create the ensembles & save them
    create_ensembles( amoptd )
    ample_util.saveAmoptd(amoptd)
