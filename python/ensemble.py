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
import logging
import os
import sys

# our imports
import ample_ensemble
import printTable
import run_spicker

def create_ensembles( amoptd ):
    """Create the ensembles using the values in the amoptd dictionary"""

    #---------------------------------------
    # Generating  ensembles
    #---------------------------------------
    logger = logging.getLogger()
    logger.info('----- Clustering models --------')
    
    # Spicker Alternative for clustering then MAX
    amoptd['spicker_rundir'] = os.path.join( amoptd['work_dir'], 'spicker_run')
    spickerer = run_spicker.SpickerCluster( amoptd )
    spickerer.run_spicker()

    logger.info( spickerer.results_summary() )

    amoptd['spicker_results'] = spickerer.results
    
    # Generate list of clusters to sample
    # Done like this so it's easy to run specific clusters
    cluster_nums = [ i for i in range(1,amoptd['num_clusters']+1)]
    
    logger.info('Clustering Done. Using the following clusters: {0}'.format(cluster_nums) )
    
    amoptd['ensemble_results'] = [] 
    for cluster in cluster_nums:
        
        logger.info('----- Truncating models for cluster {0} --------'.format(cluster)   )
        
        # Get list of models from spicker results
        cluster_models = spickerer.results[ cluster-1 ].pdb_list
        
        if len( cluster_models ) < 2:
            msg = "Could not create ensemble for cluster {0} as less than 2 models!".format( cluster )
            logger.critical(msg)
            raise RuntimeError, msg
        
        ensembler = ample_ensemble.Ensembler()
        ensembler.maxcluster_exe = amoptd['maxcluster_exe']
        ensembler.theseus_exe =  amoptd['theseus_exe']
        ensembler.max_ensemble_models = amoptd['max_ensemble_models']
        ensembler.generate_ensembles( cluster_models=cluster_models,
                                      root_dir=amoptd['work_dir'],
                                      ensemble_id=cluster,
                                      percent=amoptd['percent'] )
        # Add results to amopt
        amoptd['ensemble_results'].append( ensembler.ensembles )
        
        # List of pdbs
        #final_ensembles = ensembler.pdb_list()
        
        # Prune down to the top model        
        if amoptd['top_model_only']:
            raise RuntimeError, "Need to fix top_model_only"
            final_ensembles = truncateedit_MAX.One_model_only( final_ensembles, amoptd['work_dir'] )

        # Add to the total list
        #ensembles.append( final_ensembles )
        
        logger.info("Truncating Done for cluster {0}".format( cluster ) )
        logger.info('Created {0} ensembles'.format( len(  amoptd['ensemble_results'][-1]  ) ) )
        
    # Write out pickle file
    f = open( amoptd['results_path'], 'w' )
    cPickle.dump( amoptd, f )
    f.close()
    logging.info("Saved results as file: {0}\n".format( amoptd['results_path'] ) )
        
    return

def ensemble_summary( amoptd ):
    """Print a summary of the ensembling process"""
    
    tableFormat = printTable.Table()
    rstr=""
    
    if amoptd.has_key('spicker_results'):
        rstr = "---- Clustering Results ----\n\n"
        tdata = [ ( "Cluster", "Num Models", "Centroid" ) ]
        for i, r in enumerate( amoptd['spicker_results'] ):
            tdata.append( ( i+1, r.cluster_size, r.cluster_centroid ) )
        rstr += tableFormat.pprint_table( tdata )
    
    if amoptd.has_key('ensemble_results') and len( amoptd['ensemble_results'] ):
        rstr += "\n---- Truncation Results ----\n\n"
        # Reconstruct trunction data
        for i, clusterEnsemble in enumerate( amoptd['ensemble_results'] ):
            
            rstr += "---- Cluster {0} ----\n".format( i+1 )
            nensembles = 0
            truncation_thresholds = {}
            
            for ensemble in clusterEnsemble:
                
                nensembles += 1
                if ensemble.truncation_threshold not in truncation_thresholds.keys():
                    truncation_thresholds[ ensemble.truncation_threshold ] = ( ensemble.num_residues,
                                                                               { ensemble.radius_threshold : ensemble.num_models } )
                else:
                    # Already got this truncation level so add radius
                    if ensemble.radius_threshold not in truncation_thresholds[ ensemble.truncation_threshold ][1].keys():
                        truncation_thresholds[ ensemble.truncation_threshold ][1][ ensemble.radius_threshold ] = ensemble.num_models
            
            
            rstr += "\n"
            tdata = [ ( "Truncation Threshold", "Num Residues", "Radius Treshold", "Num Models" ) ]
            for threshold in sorted( truncation_thresholds.keys() ):
                nresidues = truncation_thresholds[ threshold ][0] 
                for i, radius in enumerate( sorted( truncation_thresholds[ threshold ][1].keys() ) ):
                    nmodels = truncation_thresholds[ threshold ][1][radius ]
                    if i == 0:
                        tdata.append( ( threshold, nresidues, radius, nmodels  ) )
                    else:
                        tdata.append( ( "", "", radius, nmodels  ) )
                    
            rstr += tableFormat.pprint_table( tdata )
            
            rstr += "\nGenerated {0} ensembles\n\n".format( nensembles )
                
    return rstr

            
if __name__ == "__main__":
    
    # This runs the ensembling starting from a pickled file containing an amopt dictionary.
    # - used when submitting the modelling jobs to a cluster
    
    if len(sys.argv) != 2 or not os.path.isfile( sys.argv[1] ):
        print "ensemble script requires the path to a pickled amopt dictionary!"
        sys.exit(1)
    
    # Get the amopt dictionary
    fpath = sys.argv[1]
    f = open( fpath, "r" )
    amoptd = cPickle.load( f )
    f.close()
    
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
    
    # Create the ensembles
    create_ensembles( amoptd )
