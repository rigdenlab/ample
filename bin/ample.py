#!/usr/bin/env ccp4-python


"""
Version 0.1
no speed buffs, no extra functionality, no optimisation
Basic Pipeline
for multi core desktops, not clusters

Issues:
1) if this script is killed, the rosetta will continue
2) Theseus can hang, this is killed after a timeout
3) is clustered by all atom RMSD to make fine clusters (option to cluster by CA only i included)
4) ASU content is the number of search models placed by MrBUMP. -- need to set this manually


================================

so...
Only data passed between different stages are the pdb files, so each stage could just start from a directory of files

"""

import os
import sys

# Test for environment variables
if not "CCP4" in sorted(os.environ.keys()):
    raise RuntimeError('CCP4 not found')

# Add the ample python folder to the PYTHONPATH
sys.path.append(os.path.join(os.environ["CCP4"], "share", "ample", "python"))
#root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
#sys.path.append( os.path.join( root, "python" ) )

# python imports
import argparse
import clusterize
import cPickle
import glob
import logging
import time

# Our imports
import add_sidechains_SCWRL
import ample_options
import ample_util
import benchmark
import ensemble
import fasta_parser
import mrbump_ensemble
import mtz_util
import nmr
import rosetta_model
import mrbump_results
import version

def main():
    """Main AMPLE routine.
    
    We require this as the multiprocessing module (only on **!!*%$$!! Windoze) requires that the main module
    can be imported. We there need ample to be a python script that can be imported, hence the main routine with
    its calling protected by the if __name__=="__main__":...
    """ 
    # get command line options
    parser = argparse.ArgumentParser( prog="AMPLE", description='Structure solution by abinitio modelling', prefix_chars="-")
    
    parser.add_argument('-alignment_file', type=str, nargs=1,
                       help='Alignment between homolog and target fastas')
    
    parser.add_argument('-all_atom', metavar='True/False', type=str, nargs=1,
                       help="Do all-atom Rosetta modelling (adds \"-return_full_atom true\" to rosetta arguments")
    
    parser.add_argument('-arpwarp_cycles', type=int, nargs=1,
                       help='The number of ArpWarp cycles to run')
    
    parser.add_argument('-ASU', type=int, nargs=1,
                       help='Manually specify the number of molecules in the asymmetric unit - sets the NMASu MRBUMP flag')
    
    parser.add_argument('-blast_dir', type=str, nargs=1,
                       help='Directory where ncbi blast is installed (binaries in expected in bin subdirectory)')
    
    parser.add_argument('-buccaneer_cycles', type=int, nargs=1,
                       help='The number of Bucanner rebuilding cycles to run')
    
    parser.add_argument('-cluster_method', type=str, nargs=1,
                       help='How to cluster the models for ensembling spicker/???')
    
    parser.add_argument('-ccp4_jobid', type=int, nargs=1,
                       help='Set the CCP4 job id - only needed when running from the CCP4 GUI')
    
    parser.add_argument('-constraints_file', type=str, nargs=1,
                       help='Residue constraints for ab initio modelling')
    
    parser.add_argument('-debug', metavar='True/False', type=str, nargs=1,
                       help='Run in debug mode (CURRENTLY UNUSED)')
    
    parser.add_argument('-domain_all_chains_pdb', type=str, nargs=1,
                       help='Fixed input to mr bump')
    
    parser.add_argument('-domain_termini_distance', type=str, nargs=1,
                       help='distance between termini for insert domains')
    
    parser.add_argument('-early_terminate', metavar='True/False', type=str, nargs=1,
                         help='Stop the run as soon as a success has been found.')
    
    parser.add_argument('-ensembler', type=str, nargs=1,
                       help='Use the Phenix ensembler')
    
    parser.add_argument('-ensembles_dir', type=str, nargs=1,
                       help='Path to directory containing existing ensembles')
    
    parser.add_argument('-fasta', type=str, nargs=1, required=True,
                       help='protein fasta file. (required)')
    
    parser.add_argument('-mr_sequence', type=str, nargs=1,
                       help="sequence file for crystal content (if different from what's given by -fasta)")
    
    parser.add_argument('-F', metavar='flag for F', type=str, nargs=1,
                       help='Flag for F column in the MTZ file')
    
    parser.add_argument('-frags_3mers', metavar='frags_3mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    parser.add_argument('-frags_9mers', metavar='frags_9mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    parser.add_argument('-FREE', metavar='flag for FREE', type=str, nargs=1,
                       help='Flag for FREE column in the MTZ file')
    
    parser.add_argument('-improve_template', metavar='improve_template', type=str, nargs=1,
                       help='Path to a template to improve - NMR, homolog' )
    
    parser.add_argument('-LGA', metavar='path_to_LGA dir', type=str, nargs=1,
                       help='pathway to LGA folder (not the exe) will use the \'lga\' executable. UNUSED')
    
    parser.add_argument('-make_frags', metavar='True/False', type=str, nargs=1,
                       help='set True to generate Rosetta 3mers and 9mers locally, False to import fragments')
    
    parser.add_argument('-make_models', metavar='True/False', type=str, nargs=1,
                       help='run rosetta modeling, set to False to import pre-made models (required if making models locally default True)')
    
    parser.add_argument('-maxcluster_exe', type=str, nargs=1,
                       help='Path to Maxcluster executable')
    
    parser.add_argument('-max_ensemble_models', type=str, nargs=1,
                       help='Maximum number of models permitted in an ensemble')
    
    parser.add_argument('-missing_domain', metavar='True/False', type=str, nargs=1,
                       help='Modelling a missing domain - requires domain_all_chains_pdb argument')
    
    parser.add_argument('-models', metavar='models', type=str, nargs=1,
                       help='Path to a folder of PDB decoys, or a tarred and gzipped/bziped, or zipped collection of decoys')
    
    parser.add_argument('-mr_keys', metavar='-mr_keys', type=str, nargs=1,
                       help='Additional keywords for MRBUMP - are passed through without editing')
    
    parser.add_argument('-name', metavar='job_name', type=str, nargs=1,
                       help='4-letter identifier for job [ampl]')
    
    parser.add_argument('-native_pdb', metavar='native_pdb', type=str, nargs=1,
                       help='Path to the crystal structure PDB for benchmarking.')
    
    parser.add_argument('-nmodels', metavar='number of models', type=int, nargs=1,
                       help='number of models to make (default: 1000)')
    
    parser.add_argument('-nr', metavar='nr', type=str, nargs=1,
                       help='Path to the NR non-redundant sequence database')
    
    parser.add_argument('-old_shelx', metavar='old_shelx', type=str, nargs=1,
                       help='old_shelx')
    
    parser.add_argument('-NMR_model_in', metavar='NMR_model_in', type=str, nargs=1,
                       help='use nmr input')
    
    parser.add_argument('-NMR_process', metavar='NMR_process', type=str, nargs=1,
                       help='number of times to process the models')
    
    parser.add_argument('-NMR_remodel_fasta', metavar='-NMR_remodel_fasta', type=str, nargs=1,
                       help='fasta_for_remodeling')
    
    parser.add_argument('-NMR_Truncate_only', metavar='True/False', type=str, nargs=1,
                       help='do no remodelling only truncate the NMR')
    
    parser.add_argument('-nproc', metavar='Number of Processors', type=int, nargs=1,
                       help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors." + \
                        "For cluster submission, this should be the number of processors on a node.")
    
    parser.add_argument('-num_clusters', type=int, nargs=1,
                       help='The number of Spicker clusters of the original decoys that will be sampled [1]')
    
    parser.add_argument('-output_pdb', type=str, nargs=1,
                       help='Name of the final result pdb to output [ample_output.pdb]')
    
    parser.add_argument('-percent', metavar='percent_truncation', type=str, nargs=1,
                       help='percent interval for truncation')
    
    parser.add_argument('-psipred_ss2', metavar='psipred file', type=str, nargs=1,
                       help='Psipred secondary structure prediction file')
    
    parser.add_argument('-phaser_kill', metavar='phaser_kill', type=int, nargs=1,
                       help='Time in minutes after which phaser will be killed (0 to leave running)')
    
    parser.add_argument('-phenix_exe', metavar='phenix_exe', type=str, nargs=1,
                       help='Path to Phenix executable')
    
    parser.add_argument('-quark_models', type=str, nargs=1,
                       help='A file containing quark models (either a tar.gz file with a pdb containing the models, or the pdb file itself)')
    
    parser.add_argument('-rg_reweight', metavar='radius of gyration reweight', type=float, nargs=1,
                       help='Set the Rosetta -rg_reweight flag to specify the radius of gyration reweight.')
    
    parser.add_argument('-rosetta_AbinitioRelax', metavar='rosetta_AbinitioRelax', type=str, nargs=1,
                       help='Path to Rosetta AbinitioRelax executable')
    
    parser.add_argument('-ROSETTA_cluster', metavar='path to Rosettas cluster', type=str, nargs=1,
                       help='location of rosetta cluster')
    
    parser.add_argument('-rosetta_db', metavar='rosetta_db', type=str, nargs=1,
                       help='Path to the Rosetta database directory')
    
    parser.add_argument('-rosetta_dir', metavar='rosetta_dir', type=str, nargs=1,
                       help='The Rosetta install directory')
    
    parser.add_argument('-rosetta_fragments_exe', metavar='rosetta_fragments_exe', type=str, nargs=1,
                       help='Location of the Rosetta make_fragments.pl script')
    
    parser.add_argument('-rosetta_version', metavar='rosetta_version', type=float, nargs=1,
                       help='The version number of Rosetta')
    
    parser.add_argument('-run_dir', metavar='run_directory', type=str, nargs=1,
                       help='Directory where the AMPLE run directory will be created [current dir]')
    
    parser.add_argument('-scwrl_exe', metavar='path to scwrl', type=str, nargs=1,
                       help='Path to Scwrl4 executable')
    
    parser.add_argument('-shelx_cycles', type=str, nargs=1,
                         help='The number of shelx cycles to run when rebuilding.')
    
    parser.add_argument('-shelxe_exe', metavar='path to shelxe executable', type=str, nargs=1,
                       help='Path to the shelxe executable')
    
    parser.add_argument('-shelxe_rebuild', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with buccaneer and arpwarp')
    
    parser.add_argument('-SIGF', metavar='SIGF', type=str, nargs=1,
                       help='Flag for SIGF column in the MTZ file')
    
    parser.add_argument('-spicker_exe', type=str, nargs=1,
                       help='Path to spicker executable')
    
    parser.add_argument('-split_mr', metavar='True/False', type=str, nargs=1,
                       help='Split MRBUMP Molecular Replacement jobs (phaser, molrep etc) into separate jobs')

    parser.add_argument('-submit_array', metavar='True/False', type=str, nargs=1,
                       help='Submit SGE jobs as array jobs')
    
    parser.add_argument('-submit_cluster', metavar='True/False', type=str, nargs=1,
                       help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')
    
    parser.add_argument('-submit_qtype', type=str, nargs=1,
                       help='cluster submission queue type - currently support SGE and LSF')
    
    parser.add_argument('-theseus_exe', metavar='Theseus exe (required)', type=str, nargs=1,
                       help='Path to theseus executable')
    
    parser.add_argument('-top_model_only', metavar='True/False', type=str, nargs=1,
                       help='Only process the top model in each ensemble')
    
    parser.add_argument('-transmembrane', metavar='transmembrane', type=str, nargs=1,
                       help='Do Rosetta modelling for transmembrane proteins')
    
    parser.add_argument('-transmembrane_octopusfile', metavar='transmembrane_octopusfile', type=str, nargs=1,
                       help='Octopus transmembrane topology predicition file')
    
    parser.add_argument('-transmembrane_spanfile', metavar='transmembrane_spanfile', type=str, nargs=1,
                       help='Span file for modelling transmembrane proteins')
    
    parser.add_argument('-transmembrane_lipofile', metavar='transmembrane_lipofile', type=str, nargs=1,
                       help='Lips4 file for modelling transmembrane proteins')

    parser.add_argument('-truncation_method', type=str, nargs=1,
                       help='How to truncate the models for ensembling percent|thresh')
    
    parser.add_argument('-truncation_pruning', type=str, nargs=1,
                       help='Whether to remove isolated residues none|single')
    
    parser.add_argument('-use_homs', metavar='True/False', type=str, nargs=1,
                       help='True =use nhomologs, False= dont use them ')
    
    parser.add_argument('-use_arpwarp', metavar='True/False', type=str, nargs=1,
                       help='True to use arpwarp to rebuild.')
    
    parser.add_argument('-use_buccaneer', metavar='True/False', type=str, nargs=1,
                       help='True to use Buccaneer')
    
    parser.add_argument('-use_scwrl', metavar='True/False', type=str, nargs=1,
                       help='Remodel sidechains of the decoy models using Scwrl4')
    
    parser.add_argument('-use_shelxe', metavar='True/False', type=str, nargs=1,
                       help='True to use shelxe')
    
    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__) )
    
    # Mutually exclusive options
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-mtz', metavar='MTZ in', type=str, nargs=1,
                       help='The MTZ file with the reflection data.')
    group.add_argument('-sf_cif', metavar='sf_cif', type=str, nargs=1,
                       help='Path to a structure factor CIF file (instead of MTZ file)')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-devel_mode', metavar='devel_mode', type=str, nargs=1,
                       help='Preset options to run in development mode - takes longer')
    group.add_argument('-quick_mode', metavar='quick_mode', type=str, nargs=1,
                       help='Preset options to run quickly, but less thoroughly')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-molrep_only', metavar='molrep_only', type=str, nargs=1,
                       help='Only use Molrep for Molecular Replacement step in MRBUMP')
    group.add_argument('-phaser_only', metavar='phaser_only', type=str, nargs=1,
                       help='Only use Phaser for Molecular Replacement step in MRBUMP')
    
    # Create amopt object to hold options
    amopt = ample_options.AmpleOptions( )
    
    # Remember the original arguments
    orig_argv = " ".join( sys.argv )
    
    # MRkeys hack - get MRBUMP keywords direct as there can be multiple arguments
    # to each one
    MRkeys = []
    #print sys.argv
    keycount = 0
    toRemove = [] # We remove all the mrkeywords
    while keycount < len(sys.argv):
        #print sys.argv[keycount] ,  keycount
        if sys.argv[keycount] == "-mr_keys":
            toRemove.append(keycount)
            keycount += 1
            tmp = []
            while not sys.argv[keycount].startswith("-"):
                #MRkeys.append( sys.argv[keycount] )
                tmp.append( sys.argv[keycount] )
                toRemove.append( keycount )
                keycount += 1
                if keycount == len(sys.argv):
                    break
            # Got last of list so add to MRkeys
            MRkeys.append( tmp )
            continue
        keycount += 1
    
    # Need to decrement the index of the items to remove
    # as we remove them
    for count, i in enumerate(toRemove):
        del sys.argv[i-count]
    ## End MRKeys hack..
    
    # convert args to dictionary
    args = parser.parse_args()
    
    # Now put them in the amopt object - this also sets/checks any defaults
    amopt.populate( args )
    
    # Now set MRKeys - might already have pre-set values so check if None or a list
    if isinstance( amopt.d['mr_keys'],list ):
        amopt.d['mr_keys'] += MRkeys
    else:
        amopt.d['mr_keys'] = MRkeys
    
    # Make a work directory and go there - this way all output goes into this directory
    if not os.path.exists( amopt.d['run_dir'] ):
        print 'Cannot find run directory: {0}'.format( amopt.d['run_dir'] )
        sys.exit()
    
    print 'Making a Run Directory: checking for previous runs\n'
    amopt.d['work_dir'] = ample_util.make_workdir( amopt.d['run_dir'], ccp4_jobid=amopt.d['ccp4_jobid'] )
    #amopt.d['work_dir'] = ample_util.make_workdir( amopt.d['run_dir'], rootname="ENSEMBLE_2_" )
    os.chdir( amopt.d['work_dir'] )
    
    # Set up logging
    logger = ample_util.setup_logging()
    
    # Print out Version and invocation
    logger.info( """AMPLE version: {0}\n\nInvoked with command-line:\n\n{1}""".format( version.__version__, orig_argv ) )
    
    # Path for pickling results
    amopt.d['results_path'] = os.path.join( amopt.d['work_dir'], "resultsd.pkl" )
    
    ###############################################################################
    #
    # FASTA processing
    #
    ###############################################################################
    # Check to see if mr_sequence was given and if not mr_sequence defaults to fasta
    if amopt.d['mr_sequence'] != None:
        if not ( os.path.exists( str(amopt.d['mr_sequence']) )):
            msg = 'Cannot find mr sequence file: {0}'.format( amopt.d['mr_sequence'] )
            logger.critical(msg)
            sys.exit(1)
    else:
        amopt.d['mr_sequence']=amopt.d['fasta']
    
    # Check we can find the input fasta
    if not ( os.path.exists( str(amopt.d['fasta']) ) or os.path.exists( str( amopt.d['NMR_remodel_fasta'] ) ) ):
        msg = 'Cannot find fasta file: {0}'.format( amopt.d['fasta'] )
        logger.critical(msg)
        sys.exit(1)
    
    # Reformat to what we need
    logger.debug('Parsing FASTA file')
    outfasta = os.path.join( amopt.d['work_dir'], amopt.d['name'] + '_.fasta')
    fp = fasta_parser.FastaParser()
    fp.reformat_fasta( amopt.d['fasta'], outfasta )
    amopt.d['fasta'] = outfasta
    amopt.d['fasta_length'] = fp.length
    logger.info( "Fasta is {0} amino acids long".format( amopt.d['fasta_length'] ) )
    
    # Check we have a decent length
    if amopt.d['fasta_length'] < 9:
        msg = "ERROR! Fasta is of length {0}. This is much too short!".format( amopt.d['fasta_length'] )
        logger.critical(msg)
        sys.exit(1)
    
    # Check we will be able to truncate at this level
    if ( float( amopt.d['fasta_length'] ) / 100 ) * float( amopt.d['percent'] ) < 1:
        msg = "Cannot truncate a fasta sequence of length {0} with {1} percent intervals. Please select a larger interval.".format( amopt.d['fasta_length'], amopt.d['percent'] )
        logger.critical(msg)
        sys.exit(1)
    
    #
    # Not sure if name actually required - see make_fragments.pl
    #
    if amopt.d['name'] and len(amopt.d['name']) != 4:
        msg = '-name argument is the wrong length, use 4 chars eg ABCD'
        logger.critical(msg)
        sys.exit(1)
        
    # Underscore required by rosetta make_fragments.pl
    amopt.d['name'] += '_'
    
    ###############################################################################
    #
    # MTZ file processing
    #
    ###############################################################################
    mtz_util.processReflectionFile( amopt.d )
    logger.info( "Using MTZ file: {0}".format( amopt.d['mtz'] ) )
    
    ###############################################################################
    #
    # Modelling and ensemble options
    #
    ###############################################################################
    
    # Set default name for modelling directory
    amopt.d['models_dir'] = amopt.d['work_dir'] + os.sep + "models"
    
    # Check if importing ensembles
    if amopt.d['ensembles_dir']:
        if not os.path.isdir( amopt.d['ensembles_dir'] ) or not len( glob.glob( os.path.join( amopt.d['ensembles_dir'], "*.pdb" ) ) ):
            msg = "Cannot import ensembles from the directory: {0}".format(amopt.d['ensembles_dir'])
            logger.critical(msg)
            sys.exit(1)
    
        amopt.d['import_ensembles'] = True
        logger.info("Found directory with ensemble files: {0}\n".format( amopt.d['ensembles_dir'] ) )
        amopt.d['make_frags'] = False
        amopt.d['make_models'] = False
    
    # Importing quark models
    if amopt.d['quark_models']:
        logger.info("Attempting to use quark models from file: {0}".format(amopt.d['quark_models']))
        if not os.path.isfile(amopt.d['quark_models']):
            logger.critical('Cannot find quark model file: {0}'.format(amopt.d['quark_models']))
            sys.exit(1)
            
        dirname,fname = os.path.split(amopt.d['quark_models'])
        # First see if we need to extract the files from the archive
        if fname.endswith(".tar.gz") or fname.endswith(".tgz"):
            logger.info("Extracting quark alldecoy.pdb file from archive: {0}".format(amopt.d['quark_models']))
            # Extract the alldecoy.pdb file from the archive
            fname=ample_util.extractFile(amopt.d['quark_models'],fileName='alldecoy.pdb',directory=amopt.d['work_dir'])
            if not fname:
                logger.critical('Cannot extract alldecoy.pdb from archive: {0}'.format(amopt.d['quark_models']))
                sys.exit(1)            
            amopt.d['quark_models']=os.path.join(amopt.d['work_dir'],fname) # Construct new path to file
    
        # Need to create the quark pdbs from the monolithic quark pdb
        logger.info("Creating directory of quark models from: {0}".format(amopt.d['quark_models']))
        quark_models=os.path.join(amopt.d['work_dir'],'quark_models')
        ample_util.splitQuark(amopt.d['quark_models'],directory=quark_models)
        amopt.d['models_dir']=quark_models
    elif amopt.d['models']:
        amopt.d['models_dir']=ample_util.extractModels(amopt.d['models'],amopt.d['models_dir'])
        amopt.d['import_models'] = True
        amopt.d['make_frags'] = False
        amopt.d['make_models'] = False
    
    # Check import flags
    if amopt.d['import_ensembles'] and (amopt.d['import_models']):
            msg = "Cannot import both models and ensembles/clusters!"
            logger.critical(msg)
            sys.exit(1)
    
    # NMR Checks
    if amopt.d['NMR_model_in']:
        if not os.path.isfile( amopt.d['NMR_model_in'] ):
            msg = "NMR_model_in flag given, but cannot find file: {0}".format( amopt.d['NMR_model_in'] )
            logger.critical(msg)
            sys.exit(1)
    
        amopt.d['NMR_protocol'] = True
        amopt.d['make_frags'] = False
        amopt.d['make_models'] = False
    
        if not os.path.isfile( str(amopt.d['NMR_remodel_fasta']) ):
            amopt.d['NMR_remodel_fasta'] =  amopt.d['fasta']
    
    if amopt.d['make_models']:
        if not os.path.isdir(amopt.d['models_dir']): os.mkdir(amopt.d['models_dir'])
        # If the user has given both fragment files we check they are ok and unset make_frags
        if amopt.d['frags_3mers'] and amopt.d['frags_9mers']:
            if not os.path.isfile( amopt.d['frags_3mers'] ) or not os.path.isfile( amopt.d['frags_9mers'] ):
                msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format( amopt.d['frags_3mers'], amopt.d['frags_9mers'] )
                logger.critical(msg)
                sys.exit(1)
            amopt.d['make_frags'] = False
    
        if amopt.d['make_frags'] and ( amopt.d['frags_3mers'] or  amopt.d['frags_9mers'] ):
            msg = "make_frags set to true, but you have given the path to the frags_3mers or frags_9mers"
            logger.critical(msg)
            sys.exit(1)
    
        if not amopt.d['make_frags'] and not ( amopt.d['frags_3mers'] and amopt.d['frags_9mers'] ):
            msg = """*** Missing fragment files! ***
    Please supply the paths to the fragment files using the -frags_3mers and -frags_9mers flags.
    These can be generated using the Robetta server: http://robetta.bakerlab.org
    Please see the AMPLE documentation for further information."""
            logger.critical(msg)
            sys.exit(1)
    
    ###############################################################################
    #
    # Misc options
    #
    ###############################################################################
    
    # Missing domains
    if amopt.d['missing_domain']:
        logger.info('Processing missing domain\n')
        if not os.path.exists( amopt.d['domain_all_chains_pdb'] ):
            msg = 'Cannot find file domain_all_chains_pdb: {0}'.format( amopt.d['domain_all_chains_pdb'] )
            logger.critical( msg )
            sys.exit(1)
    
    # MR programs
    if amopt.d['molrep_only'] and amopt.d['phaser_only']:
            logger.critical('you say you want molrep only AND phaser only, choose one or both')
            sys.exit(1)
    
    if amopt.d['molrep_only']:
        amopt.d['mrbump_programs'] = [ 'molrep' ]
    elif amopt.d['phaser_only']:
        amopt.d['mrbump_programs'] = [ 'phaser' ]
    else:
        amopt.d['mrbump_programs'] = ['molrep', 'phaser']
    
    #
    # Benchmark Mode
    #
    if amopt.d['native_pdb']:
        if not os.path.isfile(amopt.d['native_pdb']):
            msg = "Cannot find crystal structure PDB: {0}".format(amopt.d['native_pdb'])
            logger.critical(msg)
            sys.exit(1)
        amopt.d['benchmark_mode']=True
        logger.info("*** AMPLE running in benchmark mode ***")
        amopt.d['benchmark_dir']=os.path.join(amopt.d['work_dir'],"benchmark")
        os.mkdir(amopt.d['benchmark_dir'])

    ###############################################################################
    #
    # Program defaults
    #
    #
    ###############################################################################
    
    # Model building programs
    if amopt.d['use_arpwarp']:
        if not ( os.environ.has_key('warpbin') and os.path.isfile( os.path.join(os.environ['warpbin'], "auto_tracing.sh") ) ):
            logger.warn('Cannot find arpwarp script! Disabling use of arpwarp.')
            amopt.d['use_arpwarp'] = False
        else:
            logger.info('Using arpwarp script: {0}'.format( os.path.join(os.environ['warpbin'], "auto_tracing.sh") ) )
    
    #
    #Check we can find all the required programs
    #
    # Maxcluster handled differently as we may need to download the binary
    amopt.d['maxcluster_exe'] = ample_util.find_maxcluster( amopt )
    
    #
    # SPICKER and Theseus now shipped with CCP4
    #
    if not amopt.d['spicker_exe']:
        if sys.platform.startswith("win"):
            amopt.d['spicker_exe']='spicker.exe'
        else:
            amopt.d['spicker_exe']='spicker'
    try:
        amopt.d['spicker_exe'] = ample_util.find_exe(amopt.d['spicker_exe'])
    except Exception:
        logger.critical("Cannot find spicker executable: {0}".format(amopt.d['spicker_exe']))
        sys.exit(1)    
    #
    # Ensembler
    #
    if amopt.d['ensembler']:
        logger.info('You are using Phenix ensembler')
        amopt.d['phenix_exe'] = ample_util.find_exe(amopt.d['phenix_exe'])
    else:
        if not amopt.d['theseus_exe']:
            if sys.platform.startswith("win"):
                amopt.d['theseus_exe']='theseus.exe'
            else:
                amopt.d['theseus_exe']='theseus'
        try:
            amopt.d['theseus_exe'] = ample_util.find_exe(amopt.d['theseus_exe'])
        except Exception:
            logger.critical("Cannot find theseus executable: {0}".format(amopt.d['theseus_exe']))
            sys.exit(1)
    #
    # Scwrl
    #
    if amopt.d['use_scwrl']:
        if not amopt.d['scwrl_exe']:
            if sys.platform.startswith("win"):
                amopt.d['scwrl_exe']='Scwrl4.exe'
            else:
                amopt.d['scwrl_exe']='Scwrl4'
        if True:
        #try:
            amopt.d['scwrl_exe'] = ample_util.find_exe(amopt.d['scwrl_exe'])
        #except Exception:
        #    logger.critical("Cannot find Scwrl executable: {0}".format(amopt.d['scwrl_exe']))
        #    sys.exit(1)
    
    #
    # We use shelxe by default so if we can't find it we just warn and set use_shelxe to False
    #
    if amopt.d['use_shelxe']:
        if not amopt.d['shelxe_exe']:
            if sys.platform.startswith("win"):
                amopt.d['shelxe_exe']='shelxe.exe'
            else:
                amopt.d['shelxe_exe']='shelxe'
        try:
            amopt.d['shelxe_exe'] = ample_util.find_exe(amopt.d['shelxe_exe'])
        except Exception:
            msg = """*** Cannot find shelxe executable in PATH - turning off use of SHELXE. ***
    SHELXE is recommended for the best chance of success. We recommend you install shelxe from:
    http://shelx.uni-ac.gwdg.de/SHELX/
    and install it in your PATH so that AMPLE can use it.
    """
            logger.warn( msg )
            amopt.d['use_shelxe'] = False
    #
    # If shelxe_rebuild is set we need use_shelxe to be set
    #
    if amopt.d['shelxe_rebuild'] and not amopt.d['use_shelxe']:
        msg = 'shelxe_rebuild is set but use_shelxe is False. Please make sure you have shelxe installed.'
        logger.critical(msg)
        sys.exit(1)
    
    # Create the rosetta modeller - this runs all the checks required
    if amopt.d['make_models'] or amopt.d['make_frags'] or amopt.d['NMR_protocol']:  # only need Rosetta if making models
        logger.info('Using ROSETTA so checking options')
        rosetta_modeller = rosetta_model.RosettaModel(optd=amopt.d)
    
    
    if amopt.d['make_frags']:
        if amopt.d['use_homs']:
            logger.info('Making fragments (including homologues)')
        else:
            logger.info('Making fragments EXCLUDING HOMOLOGUES')
    else:
        logger.info('NOT making Fragments')
    
    if amopt.d['make_models']:
        logger.info('\nMaking Rosetta Models')
    else:
        logger.info('NOT making Rosetta Models')
        
        # Print out what is being done
    if amopt.d['use_buccaneer']:
        logger.info('Rebuilding in Bucaneer')
    else:
        logger.info('Not rebuilding in Bucaneer')
    
    if amopt.d['use_arpwarp']:
        logger.info('Rebuilding in ARP/wARP')
    else:
        logger.info('Not rebuilding in ARP/wARP')
    
    # cluster queueing
    if amopt.d['submit_cluster'] and not amopt.d['submit_qtype']:
        msg = 'Must use -submit_qtype argument to specify queueing system (e.g. QSUB, LSF ) if submitting to a cluster.'
        logger.critical(msg)
        sys.exit(1)
    
    logger.info('All needed programs are found, continuing Run')
    
    # ample_log is the 'official' output file
    ample_log = open(os.path.join(amopt.d['work_dir'],'AMPLE.log'), "w")
    ample_log.write(ample_util.header)
    ample_log.flush()
    
    # params used
    with open(os.path.join( amopt.d['work_dir'], 'params_used.txt' ), "w") as f:
        param_str = amopt.prettify_parameters()
        f.write(param_str)
    # Echo to log too
    logger.debug(param_str)
    
    #######################################################
    #
    # SCRIPT PROPER STARTS HERE
    #
    ######################################################
    time_start = time.time()
    
    # Do The Modelling
    
    # Make Rosetta fragments
    if amopt.d['make_frags']:
        rosetta_modeller.generate_fragments( submit_cluster=amopt.d['submit_cluster'],
                                             submit_qtype=amopt.d['submit_qtype'],
                                             nproc=amopt.d['nproc'] )
        amopt.d['frags_3mers'] = rosetta_modeller.frags_3mers
        amopt.d['frags_9mers'] = rosetta_modeller.frags_9mers
    
    # break here for NMR (frags needed but not modelling
    # if NMR process models first
    if amopt.d['NMR_protocol']:
        nmr.doNMR(amopt, rosetta_modeller, logger)
    # return from nmr with models already made
    elif amopt.d['make_models']:
    # Make the models
    
        logger.info('----- making Rosetta models--------')
        logger.info('making ' + str(amopt.d['nmodels']) + ' models...')
    
        # If we are running with cluster support submit all modelling jobs to the cluster queue
        if amopt.d['submit_cluster']:
            # Invoke the cluster run class
            cluster_run = clusterize.ClusterRun()
            cluster_run.modelOnCluster(rosetta_modeller, amopt.d)
        else:
            # Run locally
            amopt.d['models_dir'] = rosetta_modeller.doModelling()
    
        ##End IF amopt.d['submit_cluster']
    
        msg = 'Modelling complete - models stored in:\n   ' + amopt.d['models_dir'] + '\n'
    elif amopt.d['import_models']:
        msg = 'Importing models from directory:\n   ' + amopt.d['models_dir'] + '\n'
        ample_log.write(msg)
        logger.info(msg)
        if not ample_util.check_pdbs(amopt.d['models_dir']):
            msg = "Cannot import models from the directory: {0}".format(amopt.d['models_dir'])
            logger.critical(msg)
            sys.exit(1)
        if amopt.d['use_scwrl']:
            msg = "Processing sidechains of imported models from {0} with Scwl\n".format( amopt.d['models_dir'] )
            models_dir_scwrl = os.path.join(amopt.d['work_dir'],os.path.basename(amopt.d['models_dir'])+"_scwrl")
            if os.path.isdir( models_dir_scwrl ):
                emsg = "Scwrl models directory {0} already exists-please move it aside".format( models_dir_scwrl )
                logger.critical( emsg )
                raise RuntimeError, emsg
            os.mkdir( models_dir_scwrl )
            msg += "Scwrl-processed models will be placed in directory: {0}".format( models_dir_scwrl )
            msg += "Running Scwrl..."
            ample_log.write(msg)
            logger.info( msg )
            scwrl = add_sidechains_SCWRL.Scwrl( scwrlExe=amopt.d['scwrl_exe'] )
            scwrl.processDirectory( inDirectory=amopt.d['models_dir'], outDirectory=models_dir_scwrl )
            amopt.d['models_dir'] = models_dir_scwrl
            logger.info( "Finished processing models with Scwrl" )
    
    #---------------------------------------
    # Do the clustering
    #---------------------------------------
    
    ## Save results
    #f = open( amopt.d['results_path'], 'w' )
    #cPickle.dump( amopt.d, f )
    #f.close()
    #logging.info("Saved results as file: {0}\n".format( amopt.d['results_path'] ) )
    #
    ensembles = [] # List of ensembles - 1 per cluster
    if amopt.d['import_ensembles']:
        # Importing pre-made ensembles
    
        # Set list of ensembles to the one we are importing
        msg = '\nImporting ensembles from directory:\n   ' + amopt.d['ensembles_dir'] + '\n\n'
        ample_log.write(msg)
        logger.info( msg )
        ensembles =  glob.glob( os.path.join(amopt.d['ensembles_dir'], '*.pdb') )
    else:
        # Check we have some models to work with
        if not glob.glob(os.path.join(amopt.d['models_dir'],"*.pdb")):
            msg="ERORR! Cannot find any pdb files in: {0}".format(amopt.d['models_dir'])
            logger.critical(msg)
            ample_util.saveAmoptd(amopt.d)
            sys.exit(1)
        
        if amopt.d['submit_cluster']:
    
            # Pickle dictionary so it can be opened by the job to get the parameters
            with open( amopt.d['results_path'], 'w' ) as f:
                cPickle.dump( amopt.d, f )
    
            mrBuild = clusterize.ClusterRun()
            mrBuild.ensembleOnCluster( amopt.d )
            mrBuild.monitorQueue()
    
            # queue finished so unpickle results
            with open( amopt.d['results_path'], "r" ) as f: amopt.d = cPickle.load( f )
        else:
            ensemble.create_ensembles( amopt.d )
    
        # Check we have something to work with
        if not amopt.d.has_key('ensembles') or not len( amopt.d['ensembles'] ):
            msg = "Could not load any ensembles after running create_ensembles!"
            logger.critical( msg )
            sys.exit(1)
            
        ensembles=amopt.d['ensembles']
        ensemble_summary = ensemble.ensemble_summary(amopt.d['ensembles_data'])
        ample_log.write(ensemble_summary)
        ample_log.flush()
        logger.info(ensemble_summary)
    #
    # Bail here if we didn't create anything
    #
    if not len(ensembles):
        logger.critical("### AMPLE FAILED TO GENERATE ANY ENSEMBLES! ###\nExiting...")
        ample_util.saveAmoptd(amopt.d)
        sys.exit(1)
    
    # MRBUMP analysis of the ensembles
    logger.info('----- Running MRBUMP on ensembles--------\n\n')
    if len(ensembles) < 1:
        msg = "ERROR! Cannot run MRBUMP as there are no ensembles!"
        logger.critical( msg )
        sys.exit()
    
    bump_dir = os.path.join(amopt.d['work_dir'], 'MRBUMP')
    if not os.path.exists(bump_dir):
        os.mkdir(bump_dir)
    os.chdir(bump_dir)
    amopt.d['mrbump_dir'] = bump_dir
    amopt.d['mrbump_results'] = []
    logger.info("Running MRBUMP jobs in directory: {0}".format( bump_dir ) )
    
    # Create job scripts
    logger.info("Generating MRBUMP runscripts")
    job_scripts = mrbump_ensemble.generate_jobscripts( ensembles, amopt.d )
    #continue

    if amopt.d['submit_cluster']:
        mrbump_ensemble.mrbump_ensemble_cluster( job_scripts, amopt.d )
    else:
        mrbump_ensemble.mrbump_ensemble_local( job_scripts, amopt.d )

    # Collect the MRBUMP results
    results_summary = mrbump_results.ResultsSummary()
    results_summary.extractResults(bump_dir)
    amopt.d['mrbump_results'] = results_summary.results
    ample_util.saveAmoptd(amopt.d)
    
    # Timing data
    time_stop = time.time()
    elapsed_time = time_stop - time_start
    run_in_min = elapsed_time / 60
    run_in_hours = run_in_min / 60
    msg = '\nMR and shelx DONE\n\n ALL DONE  (in ' + str(run_in_hours) + ' hours) \n----------------------------------------\n'
    logging.info(msg)
    ample_log.write(msg)
    ample_log.flush()
    
    # Benchmark mode
    if amopt.d['benchmark_mode']:
        benchmark.analyse(amopt.d)
        ample_util.saveAmoptd(amopt.d)

    # Now print out the final summary
    summary = mrbump_results.finalSummary(amopt.d)
    logger.info( summary )
    ample_log.write(summary)
    ample_log.flush()
    ample_log.close()

    sys.exit(0)
    #------------------------------------
    # END
    #----------------------------------

if __name__=="__main__":
    main()

