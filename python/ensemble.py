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

# our imports
import ample_ensemble
import ample_util
import printTable
import run_spicker


def spicker_cluster( amoptd ):

    logger = logging.getLogger()
    logger.info('----- Clustering models --------')

    # Spicker Alternative for clustering
    amoptd['spicker_rundir'] = os.path.join( amoptd['work_dir'], 'spicker_run')
    spickerer = run_spicker.SpickerCluster( amoptd )
    spickerer.run_spicker()
    logger.info( spickerer.results_summary() )
    amoptd['spicker_results'] = spickerer.results

    return

def ensemble_models( cluster_models, amoptd, ensemble_id='X' ):

    logger = logging.getLogger()

    logger.info('----- Truncating models for cluster {0} --------'.format( ensemble_id ) )

    if len( cluster_models ) < 2:
        msg = "Could not create ensemble for cluster {0} as less than 2 models!".format( ensemble_id )
        logger.critical(msg)
        #raise RuntimeError, msg
        return None

    ensembler = ample_ensemble.Ensembler()
    ensembler.maxcluster_exe = amoptd['maxcluster_exe']
    ensembler.theseus_exe =  amoptd['theseus_exe']
    ensembler.max_ensemble_models = amoptd['max_ensemble_models']
    ensembler.generate_ensembles( cluster_models=cluster_models,
                                  root_dir=amoptd['work_dir'],
                                  ensemble_id=ensemble_id,
                                  percent=amoptd['percent'],
                                  mode=amoptd['ensemble_mode']
                                   )

    return ensembler.ensembles

def create_ensembles( amoptd ):
    """Create the ensembles using the values in the amoptd dictionary"""

    logger = logging.getLogger()

    #---------------------------------------
    # Generating  ensembles
    #---------------------------------------

    # Cluster with Spicker
    spicker_cluster( amoptd )

    # Generate list of clusters to sample
    # Done like this so it's easy to run specific clusters
    cluster_nums = [ i for i in range(1,amoptd[ 'num_clusters' ]+1 )]

    logger.info('Clustering Done. Using the following clusters: {0}'.format( cluster_nums ) )

    amoptd['ensemble_results'] = []
    for cluster in cluster_nums:

        # Get list of models from spicker results
        ensembles = ensemble_models( amoptd['spicker_results'][ cluster-1 ].pdb_list, amoptd, ensemble_id=cluster )

        # Add results to amopt
        amoptd['ensemble_results'].append( ensembles )

        # Prune down to the top model
        if amoptd['top_model_only']:
            raise RuntimeError, "Need to fix top_model_only"
            #final_ensembles = truncateedit_MAX.One_model_only( final_ensembles, amoptd['work_dir'] )

        logger.info("Truncating Done for cluster {0}".format( cluster ) )
        if amoptd['ensemble_results'][-1] and len(amoptd['ensemble_results'][-1]):
            logger.info('Created {0} ensembles'.format( len(  amoptd['ensemble_results'][-1] ) ) )

    return

def import_cluster( amoptd ):

    logger = logging.getLogger()

    cluster_src = amoptd['import_cluster']
    if not cluster_src:
        msg = "import_cluster no cluster set!"
        logger.critical( msg )
        raise RuntimeError, msg

    # Either read from list of files or a directory
    clusterPdbs = []
    if os.path.isdir( cluster_src ):
        logger.info( "Importing existing cluster from directory: {0}".format( cluster_src ) )
        clusterPdbs = glob.glob( cluster_src + "*.pdb" )
    elif os.path.isfile( cluster_src ):
        logger.info( "Importing existing cluster from PDBs specified in file: {0}".format( cluster_src ) )
        clusterPdbs = [ line.strip() for line in open( cluster_src ) ]
    else:
        msg = "import_cluster cannot import cluster from: {0}".format( cluster_src )
        logger.critical( msg )
        raise RuntimeError, msg

    if not len( clusterPdbs ):
        msg = "import_cluster cannot find any pdb files in: {0}".format( cluster_src )
        logger.critical( msg )
        raise RuntimeError, msg

    ensembles = ensemble_models( clusterPdbs, amoptd, ensemble_id='X' )

    amoptd['ensemble_results'] = [ ensembles ]

    return [ e.pdb for e in ensembles ]

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

            if not clusterEnsemble:
                rstr += "\n### COULD NOT GENERATE ANY ENSEMBLES! ###\n\n"
                continue

            nensembles = 0
            truncation_levels = {}

            side_chain_treatments = []
            for ensemble in clusterEnsemble:

                nensembles += 1

                if ensemble.side_chain_treatment not in side_chain_treatments:
                    side_chain_treatments.append( ensemble.side_chain_treatment )

                if ensemble.truncation_level not in truncation_levels.keys():
                    truncation_levels[ ensemble.truncation_level ] = ( ensemble.truncation_threshold, ensemble.num_residues,
                                                                               { ensemble.radius_threshold : ensemble.num_models } )
                else:
                    # Already got this truncation level so add radius
                    if ensemble.radius_threshold not in truncation_levels[ ensemble.truncation_level ][2].keys():
                        truncation_levels[ ensemble.truncation_level ][2][ ensemble.radius_threshold ] = ensemble.num_models


            rstr += "\n"
            tdata = [ ( "Truncation Level", "Variance Threshold (A^2)", "No. Residues", "Radius Threshold (A)", "No. Decoys" ) ]
            raw_ensemble_count=0
            for level in sorted( truncation_levels.keys() ):
                #truncation_level=len(truncation_levels)-level
                threshold = truncation_levels[ level ][0]
                nresidues = truncation_levels[ level ][1]
                for i, radius in enumerate( sorted( truncation_levels[ level ][2].keys() ) ):
                    raw_ensemble_count+=1
                    nmodels = truncation_levels[ level ][2][radius ]
                    if i == 0:
                        tdata.append( ( level, threshold, nresidues, radius, nmodels  ) )
                    else:
                        tdata.append( ( "", "", "", radius, nmodels  ) )

            rstr += tableFormat.pprint_table( tdata )

            rstr += "\nGenerated {0} raw ensembles\n\n".format( raw_ensemble_count )

            rstr += "Each raw ensemble will be subject to the following {0} sidechain treatments: {1}\n\n".format( len(side_chain_treatments),  ", ".join(side_chain_treatments) )
            rstr += "{0} ensembles in total have been generated\n\n".format( nensembles )

            assert raw_ensemble_count * len(side_chain_treatments) == nensembles, "Error generating correct number of ensembles!"

    return rstr


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
