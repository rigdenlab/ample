#!/usr/bin/env ccp4-python
"""
This is AMPLE

This script is named ample.py due to a problem with running the multiprocessing (which is used to parallelise
the running of jobs on a local machine - see python/workers.py) module under windows.
The multiprocessing module on windows requires that it can import the main module, and the import machinery
requires that any file being imported is named <foo>.py, and any changes to this would require hacking the 
multiprocessing module, so to avoid this, our script must be called ample.py
"""
import os
import platform
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
import cPickle
import glob
import logging
import shutil
import time

# Our imports
import add_sidechains_SCWRL
import ample_exit
import ample_options
import ample_sequence
import ample_util
import benchmark
import ensemble
import mrbump_ensemble
import mrbump_results
import mtz_util
import nmr
import pdb_edit
import pyrvapi_results
import rosetta_model
import version
import workers

def process_command_line():
    # get command line options
    parser = argparse.ArgumentParser(prog="AMPLE", description='Structure solution by abinitio modelling', prefix_chars="-")
    
    parser.add_argument('-alignment_file', type=str, nargs=1,
                       help='Alignment file in fasta format. For homologues the first line of each sequence must be the pdb file name')
    
    parser.add_argument('-all_atom', metavar='True/False', type=str, nargs=1,
                       help="Do all-atom Rosetta modelling (adds \"-return_full_atom true\" to rosetta arguments")
    
    parser.add_argument('-arpwarp_cycles', type=int, nargs=1,
                       help='The number of ArpWarp cycles to run')
    
    parser.add_argument('-blast_dir', type=str, nargs=1,
                       help='Directory where ncbi blast is installed (binaries in expected in bin subdirectory)')
    
    parser.add_argument('-buccaneer_cycles', type=int, nargs=1,
                       help='The number of Bucanner rebuilding cycles to run')
    
    parser.add_argument('-cluster_dir', type=str, nargs=1,
                       help='Path to directory of pre-clustered models to import')
    
    parser.add_argument('-cluster_method', type=str, nargs=1,
                       help='How to cluster the models for ensembling (spicker|fast_protein_cluster')
    
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
    
    parser.add_argument('-dry_run', metavar='True/False', type=str, nargs=1,
                         help='Check if input files and supplied options are valid.')
    
    parser.add_argument('-early_terminate', metavar='True/False', type=str, nargs=1,
                         help='Stop the run as soon as a success has been found.')
    
    parser.add_argument('-ensembles_dir', type=str, nargs=1,
                       help='Path to directory containing existing ensembles')
    
    parser.add_argument('-fasta', type=str, nargs=1, required=True,
                       help='protein fasta file. (required)')
    
    parser.add_argument('-fast_protein_cluster_exe', type=str, nargs=1,
                       help='path to fast_protein_cluster executable')
    
    parser.add_argument('-F', metavar='flag for F', type=str, nargs=1,
                       help='Flag for F column in the MTZ file')
    
    parser.add_argument('-frags_3mers', metavar='frags_3mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    parser.add_argument('-frags_9mers', metavar='frags_9mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    parser.add_argument('-FREE', metavar='flag for FREE', type=str, nargs=1,
                       help='Flag for FREE column in the MTZ file')
    
    parser.add_argument('-gesamt_exe', metavar='gesamt_exe', type=str, nargs=1,
                       help='Path to the gesamt executable')
    
    parser.add_argument('-homologs', metavar='True/False', type=str, nargs=1,
                       help='Generate ensembles from homologs models (requires -alignment_file)')
    
    parser.add_argument('-homolog_aligner', metavar='homolog_aligner', type=str, nargs=1,
                       help='Program to use for structural alignment of homologs (gesamt|mustang)')
    
    parser.add_argument('-ideal_helices', metavar='True/False', type=str, nargs=1,
                       help='Use ideal polyalanine helices to solve structure (8 helices: from 5-40 residues)')

    parser.add_argument('-improve_template', metavar='improve_template', type=str, nargs=1,
                       help='Path to a template to improve - NMR, homolog')
    
    parser.add_argument('-LGA', metavar='path_to_LGA dir', type=str, nargs=1,
                       help='pathway to LGA folder (not the exe) will use the \'lga\' executable. UNUSED')
    
    parser.add_argument('-make_frags', metavar='True/False', type=str, nargs=1,
                       help='set True to generate Rosetta 3mers and 9mers locally, False to import fragments')
    
    parser.add_argument('-make_models', metavar='True/False', type=str, nargs=1,
                       help='run rosetta modeling, set to False to import pre-made models (required if making models locally default True)')
    
    parser.add_argument('-maxcluster_exe', type=str, nargs=1,
                       help='Path to Maxcluster executable')
    
    parser.add_argument('-max_array_jobs', type=str, nargs=1,
                       help='Maximum number of array jobs to run')
    
    parser.add_argument('-max_ensemble_models', type=str, nargs=1,
                       help='Maximum number of models permitted in an ensemble')
    
    parser.add_argument('-missing_domain', metavar='True/False', type=str, nargs=1,
                       help='Modelling a missing domain - requires domain_all_chains_pdb argument')
    
    parser.add_argument('-models', metavar='models', type=str, nargs=1,
                       help='Path to a folder of PDB decoys, or a tarred and gzipped/bziped, or zipped collection of decoys')
    
    parser.add_argument('-mr_keys', metavar='-mr_keys', type=str, nargs=1,
                       help='Additional keywords for MRBUMP - are passed through without editing')

    parser.add_argument('-mr_sequence', type=str, nargs=1,
                       help="sequence file for crystal content (if different from what's given by -fasta)")

    parser.add_argument('-mustang_exe', metavar='mustang_exe', type=str, nargs=1,
                       help='Path to the mustang executable')

    parser.add_argument('-name', metavar='job_name', type=str, nargs=1,
                       help='4-letter identifier for job [ampl]')
    
    parser.add_argument('-native_pdb', metavar='native_pdb', type=str, nargs=1,
                       help='Path to the crystal structure PDB for benchmarking.')

    parser.add_argument('-nmasu', type=int, nargs=1,
                       help='Manually specify the number of molecules in the asymmetric unit - sets the NMASu MRBUMP flag')
    
    parser.add_argument('-nmodels', metavar='number of models', type=int, nargs=1,
                       help='number of models to make (default: 1000)')
    
    parser.add_argument('-nr', metavar='nr', type=str, nargs=1,
                       help='Path to the NR non-redundant sequence database')
    
    parser.add_argument('-nmr_model_in', metavar='nmr_model_in', type=str, nargs=1,
                       help='PDB with NMR models')
    
    parser.add_argument('-nmr_process', metavar='nmr_process', type=int, nargs=1,
                       help='number of times to process the NMR models')
    
    parser.add_argument('-nmr_remodel', metavar='True/False', type=str, nargs=1,
                       help='Remodel the NMR structures')
    
    parser.add_argument('-nmr_remodel_fasta', metavar='nmr_remodel_fasta', type=str, nargs=1,
                       help='The FASTA sequence to be used for remodelling the NMR ensemble if different from the default FASTA sequence')
    
    parser.add_argument('-no_gui', metavar='True/False', type=str, nargs=1,
                       help='Do not display the AMPLE gui.')
    
    parser.add_argument('-nproc', metavar='Number of Processors', type=int, nargs=1,
                       help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors." + \
                        "For cluster submission, this should be the number of processors on a node.")
    
    parser.add_argument('-num_clusters', type=int, nargs=1,
                       help='The number of Spicker clusters of the original decoys that will be sampled [1]')
    
    parser.add_argument('-output_pdb', type=str, nargs=1,
                       help='Name of the final result pdb to output [ample_output.pdb]')
    
    parser.add_argument('-purge', metavar='True/False', type=str, nargs=1,
                       help='Delete all intermediate files and failed MRBUMP results')
    
    parser.add_argument('-percent', metavar='percent_truncation', type=str, nargs=1,
                       help='percent interval for truncation')
    
    parser.add_argument('-psipred_ss2', metavar='psipred file', type=str, nargs=1,
                       help='Psipred secondary structure prediction file')
    
    parser.add_argument('-phaser_kill', metavar='phaser_kill', type=int, nargs=1,
                       help='Time in minutes after which phaser will be killed (0 to leave running)')
    
    parser.add_argument('-phaser_rms', metavar='phaser_rms', type=float, nargs=1,
                       help='rms value for phaser (default=0.1)')
    
    parser.add_argument('-phenix_exe', metavar='phenix_exe', type=str, nargs=1,
                       help='Path to Phenix executable')
    
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
    
    parser.add_argument('-shelxe_rebuild_arpwarp', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with arpwarp')
    
    parser.add_argument('-shelxe_rebuild_buccaneer', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with buccaneer')
    
    parser.add_argument('-SIGF', metavar='SIGF', type=str, nargs=1,
                       help='Flag for SIGF column in the MTZ file')
    
    parser.add_argument('-spicker_exe', type=str, nargs=1,
                       help='Path to spicker executable')
    
    parser.add_argument('-submit_array', metavar='True/False', type=str, nargs=1,
                       help='Submit SGE jobs as array jobs')
    
    parser.add_argument('-submit_cluster', metavar='True/False', type=str, nargs=1,
                       help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')
    
    parser.add_argument('-submit_qtype', type=str, nargs=1,
                       help='cluster submission queue type - currently support SGE and LSF')
    
    parser.add_argument('-submit_queue', type=str, nargs=1,
                       help='The queue to submit to on the cluster.')
    
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
                       help='How to truncate the models for ensembling percent|thresh|focussed')
    
    parser.add_argument('-truncation_pruning', type=str, nargs=1,
                       help='Whether to remove isolated residues none|single')
    
    parser.add_argument('-use_homs', metavar='True/False', type=str, nargs=1,
                       help='Select ROSETTA fragments from homologous models')
    
    parser.add_argument('-use_arpwarp', metavar='True/False', type=str, nargs=1,
                       help='True to use arpwarp to rebuild.')
    
    parser.add_argument('-use_buccaneer', metavar='True/False', type=str, nargs=1,
                       help='True to use Buccaneer')
    
    parser.add_argument('-use_scwrl', metavar='True/False', type=str, nargs=1,
                       help='Remodel sidechains of the decoy models using Scwrl4')
    
    parser.add_argument('-use_shelxe', metavar='True/False', type=str, nargs=1,
                       help='True to use shelxe')
    
    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__))
    
    parser.add_argument('-webserver_uri', type=str, nargs=1,
                       help='URI of the webserver directory - also indicates we are running as a webserver')
    
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
    amopt = ample_options.AmpleOptions()
    
    # Remember the original arguments
    orig_argv = " ".join(sys.argv)
    
    # MRkeys hack - get MRBUMP keywords direct as there can be multiple arguments
    # to each one
    MRkeys = []
    # print sys.argv
    keycount = 0
    toRemove = []  # We remove all the mrkeywords
    while keycount < len(sys.argv):
        # print sys.argv[keycount] ,  keycount
        if sys.argv[keycount] == "-mr_keys":
            toRemove.append(keycount)
            keycount += 1
            tmp = []
            while not sys.argv[keycount].startswith("-"):
                # MRkeys.append( sys.argv[keycount] )
                tmp.append(sys.argv[keycount])
                toRemove.append(keycount)
                keycount += 1
                if keycount == len(sys.argv):
                    break
            # Got last of list so add to MRkeys
            MRkeys.append(tmp)
            continue
        keycount += 1
    
    # Need to decrement the index of the items to remove
    # as we remove them
    for count, i in enumerate(toRemove):
        del sys.argv[i - count]
    # # End MRKeys hack..
    
    # convert args to dictionary
    args = parser.parse_args()
    
    # Now put them in the amopt object - this also sets/checks any defaults
    amopt.populate(args)

    # Now set MRKeys - might already have pre-set values so check if None or a list
    if isinstance(amopt.d['mr_keys'], list):
        amopt.d['mr_keys'] += MRkeys
    else:
        amopt.d['mr_keys'] = MRkeys
    
    return amopt, orig_argv

def process_options(amoptd, logger):
    
    # Path for pickling results
    amoptd['results_path'] = os.path.join(amoptd['work_dir'], "resultsd.pkl")
    
    ###############################################################################
    #
    # FASTA processing
    #
    ###############################################################################
    # Check to see if mr_sequence was given and if not mr_sequence defaults to fasta
    if amoptd['mr_sequence'] != None:
        if not (os.path.exists(str(amoptd['mr_sequence']))):
            msg = 'Cannot find mr sequence file: {0}'.format(amoptd['mr_sequence'])
            ample_exit.exit(msg)
    else:
        amoptd['mr_sequence'] = amoptd['fasta']
    
    # Check we can find the input fasta
    if not os.path.exists(str(amoptd['fasta'])):
        msg = 'Cannot find fasta file: {0}'.format(amoptd['fasta'])
        ample_exit.exit(msg)
    
    # Reformat to what we need
    logger.debug('Parsing FASTA file')
    fp = ample_sequence.Sequence(fasta=amoptd['fasta'])
    if fp.numSequences() != 1:
        msg = "ERROR! Fasta file {0} has > 1 sequence in it.".format(amoptd['fasta'])
        ample_exit.exit(msg)
    
    # Length checks
    amoptd['fasta_length'] = fp.length()
    logger.info("Fasta is {0} amino acids long".format(amoptd['fasta_length']))
    
    # Check we have a decent length
    if amoptd['fasta_length'] < 9:
        msg = "ERROR! Fasta is of length {0}. This is much too short!".format(amoptd['fasta_length'])
        ample_exit.exit(msg)
    
    # Check we will be able to truncate at this level
    if (float(amoptd['fasta_length']) / 100) * float(amoptd['percent']) < 1:
        msg = "Cannot truncate a fasta sequence of length {0} with {1} percent intervals. Please select a larger interval.".format(amoptd['fasta_length'], amoptd['percent'])
        ample_exit.exit(msg)
    
    # Fasta is ok, so write out a canonical fasta in the work directory
    outfasta = os.path.join(amoptd['work_dir'], amoptd['name'] + '_.fasta')
    fp.write_fasta(outfasta)
    amoptd['fasta'] = outfasta
    amoptd['sequence'] = fp.sequence()
    #
    # Not sure if name actually required - see make_fragments.pl
    #
    if amoptd['name'] and len(amoptd['name']) != 4:
        msg = '-name argument is the wrong length, use 4 chars eg ABCD'
        ample_exit.exit(msg)
        
    # Underscore required by rosetta make_fragments.pl
    amoptd['name'] += '_'
    
    ###############################################################################
    #
    # MTZ file processing
    #
    ###############################################################################
    try:
        mtz_util.processReflectionFile(amoptd)
    except Exception, e:
        ex_type, ex, tb = sys.exc_info()
        msg = "Error processing reflection file: {0}".format(e)
        ample_exit.exit(msg, tb)
    logger.info("Using MTZ file: {0}".format(amoptd['mtz']))
    
    ###############################################################################
    #
    # Modelling and ensemble options
    #
    ###############################################################################
    
    # Set default name for modelling directory
    amoptd['models_dir'] = os.path.join(amoptd['work_dir'], "models")
    
    # Check if importing ensembles
    if amoptd['ensembles_dir']:
        if not pdb_edit.check_pdb_directory(amoptd['ensembles_dir'], single=False):
            msg = "Cannot import ensembles from the directory: {0}".format(amoptd['ensembles_dir'])
            ample_exit.exit(msg)
        amoptd['import_ensembles'] = True
        logger.info("Found directory with ensemble files: {0}\n".format(amoptd['ensembles_dir']))
        amoptd['make_frags'] = False
        amoptd['make_models'] = False
    elif amoptd['cluster_dir']:
        if not os.path.isdir(amoptd['cluster_dir']):
            msg = "Import cluster cannot find directory: {0}".format(amoptd['cluster_dir'])
            ample_exit.exit(msg)
        if not glob.glob(os.path.join(amoptd['cluster_dir'], "*.pdb")):
            msg = "Import cluster cannot find pdbs in directory: {0}".format(amoptd['cluster_dir'])
            ample_exit.exit(msg)
        logger.info("Importing pre-clustered models from directory: {0}\n".format(amoptd['cluster_dir']))   
        amoptd['import_cluster'] = True
        amoptd['make_frags'] = False
        amoptd['make_models'] = False
    elif amoptd['ideal_helices']:
        amoptd['make_frags'] = False
        amoptd['make_models'] = False
    elif amoptd['homologs']:
        amoptd['make_frags'] = False
        amoptd['make_models'] = False
        if not (os.path.isfile(str(amoptd['alignment_file'])) or os.path.isfile(str(amoptd['mustang_exe'])) or os.path.isfile(str(amoptd['gesamt_exe']))):
            msg = "Homologs option requires an aligment file or path to a mustang or gesamt executable to be supplied\n" + \
            "Please supply the path to mustang or gesamt executable with the -mustang_exe or -gesamt_exe flag, or an alignment file in fasta format with the -alignment_file flag"
            ample_exit.exit(msg)
        if not os.path.isdir(str(amoptd['models'])):
            msg = "Homologs option requires a directory of pdb models to be supplied\n" + \
            "Please supply the models with the -models flag"
            ample_exit.exit(msg)
        amoptd['import_models'] = True
        amoptd['models_dir'] = ample_util.extract_models(amoptd['models'],
                                                         directory=amoptd['models_dir'],
                                                         sequence=None,
                                                         single=False,
                                                         allsame=False)
    elif amoptd['models']:
        amoptd['models_dir'] = ample_util.extract_models(amoptd['models'], amoptd['models_dir'])
        amoptd['import_models'] = True
        amoptd['make_frags'] = False
        amoptd['make_models'] = False
        
    # Check import flags
    if amoptd['import_ensembles'] and (amoptd['import_models']):
            msg = "Cannot import both models and ensembles/clusters!"
            ample_exit.exit(msg)
    
    # NMR Checks
    if amoptd['nmr_model_in']:
        msg = "Using nmr_model_in file: {0}".format(amoptd['nmr_model_in'])
        logger.info(msg)
        if not os.path.isfile(amoptd['nmr_model_in']):
            msg = "nmr_model_in flag given, but cannot find file: {0}".format(amoptd['nmr_model_in'])
            ample_exit.exit(msg)
        if amoptd['nmr_remodel']:
            if amoptd['nmr_remodel_fasta']:
                if not os.path.isfile(amoptd['nmr_remodel_fasta']):
                    msg = "Cannot find nmr_remodel_fasta file: {0}".format(amoptd['nmr_remodel_fasta'])
                    ample_exit.exit(msg)
            else:
                amoptd['nmr_remodel_fasta'] = amoptd['fasta']
            msg = "NMR model will be remodelled with ROSETTA using the sequence from: {0}".format(amoptd['nmr_remodel_fasta'])
            logger.info(msg)
            if not amoptd['frags_3mers'] and amoptd['frags_9mers']:
                amoptd['make_frags'] = True
                msg = "nmr_remodel - will be making our own fragment files"
                logger.info(msg)
            else:
                if not os.path.isfile(amoptd['frags_3mers']) or not os.path.isfile(amoptd['frags_9mers']):
                    msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(amoptd['frags_3mers'], amoptd['frags_9mers'])
                    ample_exit.exit(msg)
                amoptd['make_frags'] = False

        else:
            amoptd['make_models'] = False
            msg = "Running in NMR truncate only mode"
            logger.info(msg)

    elif amoptd['make_models']:
        if not os.path.isdir(amoptd['models_dir']): os.mkdir(amoptd['models_dir'])
        # If the user has given both fragment files we check they are ok and unset make_frags
        if amoptd['frags_3mers'] and amoptd['frags_9mers']:
            if not os.path.isfile(amoptd['frags_3mers']) or not os.path.isfile(amoptd['frags_9mers']):
                msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(amoptd['frags_3mers'], amoptd['frags_9mers'])
                ample_exit.exit(msg)
            amoptd['make_frags'] = False
        if amoptd['make_frags'] and (amoptd['frags_3mers'] or  amoptd['frags_9mers']):
            msg = "make_frags set to true, but you have given the path to the frags_3mers or frags_9mers"
            ample_exit.exit(msg)
    
        if not amoptd['make_frags'] and not (amoptd['frags_3mers'] and amoptd['frags_9mers']):
            msg = """*** Missing fragment files! ***
    Please supply the paths to the fragment files using the -frags_3mers and -frags_9mers flags.
    These can be generated using the Robetta server: http://robetta.bakerlab.org
    Please see the AMPLE documentation for further information."""
            ample_exit.exit(msg)
    
    ###############################################################################
    #
    # Misc options
    #
    ###############################################################################
    
    # Missing domains
    if amoptd['missing_domain']:
        logger.info('Processing missing domain\n')
        if not os.path.exists(amoptd['domain_all_chains_pdb']):
            msg = 'Cannot find file domain_all_chains_pdb: {0}'.format(amoptd['domain_all_chains_pdb'])
            ample_exit.exit(msg)
    # MR programs
    if amoptd['molrep_only']:
            amoptd['phaser_only'] = False
            #msg = 'you say you want molrep only AND phaser only, choose one or both'
            #ample_exit.exit(msg)
    
    if amoptd['molrep_only']:
        amoptd['mrbump_programs'] = [ 'molrep' ]
    elif amoptd['phaser_only']:
        amoptd['mrbump_programs'] = [ 'phaser' ]
    else:
        amoptd['mrbump_programs'] = ['molrep', 'phaser']
    #
    # Benchmark Mode
    #
    if amoptd['native_pdb']:
        if not os.path.isfile(amoptd['native_pdb']):
            msg = "Cannot find crystal structure PDB: {0}".format(amoptd['native_pdb'])
            ample_exit.exit(msg)
        amoptd['benchmark_mode'] = True
        logger.info("*** AMPLE running in benchmark mode ***")
        amoptd['benchmark_dir'] = os.path.join(amoptd['work_dir'], "benchmark")
        os.mkdir(amoptd['benchmark_dir'])

    ###############################################################################
    #
    # Program defaults
    #
    #
    ###############################################################################
    
    # Model building programs
    if amoptd['use_arpwarp']:
        if not (os.environ.has_key('warpbin') and os.path.isfile(os.path.join(os.environ['warpbin'], "auto_tracing.sh"))):
            logger.warn('Cannot find arpwarp script! Disabling use of arpwarp.')
            amoptd['use_arpwarp'] = False
        else:
            logger.info('Using arpwarp script: {0}'.format(os.path.join(os.environ['warpbin'], "auto_tracing.sh")))
    #
    # Check we can find all the required programs
    #
    # Maxcluster handled differently as we may need to download the binary
    amoptd['maxcluster_exe'] = ample_util.find_maxcluster(amoptd)
    
    #
    # Ensemble options
    #
    if amoptd['cluster_method'] == 'spicker':
        if not amoptd['spicker_exe']:
            if sys.platform.startswith("win"):
                amoptd['spicker_exe'] = 'spicker.exe'
            else:
                amoptd['spicker_exe'] = 'spicker'
        try:
            amoptd['spicker_exe'] = ample_util.find_exe(amoptd['spicker_exe'])
        except Exception:
            msg = "Cannot find spicker executable: {0}".format(amoptd['spicker_exe'])
            ample_exit.exit(msg)
    elif amoptd['cluster_method'] == 'fast_protein_cluster':
        if not amoptd['fast_protein_cluster_exe']: amoptd['fast_protein_cluster_exe'] = 'fast_protein_cluster'
        try:
            amoptd['fast_protein_cluster_exe'] = ample_util.find_exe(amoptd['fast_protein_cluster_exe'])
        except Exception:
            msg = "Cannot find fast_protein_cluster executable: {0}".format(amoptd['fast_protein_cluster_exe'])
            ample_exit.exit(msg)
    else:
        msg = "Unrecognised cluster_method: {0}".format(amoptd['cluster_method'])
        ample_exit.exit(msg)
    if not amoptd['theseus_exe']:
        if sys.platform.startswith("win"):
            amoptd['theseus_exe'] = 'theseus.exe'
        else:
            amoptd['theseus_exe'] = 'theseus'
    try:
        amoptd['theseus_exe'] = ample_util.find_exe(amoptd['theseus_exe'])
    except Exception:
        msg = "Cannot find theseus executable: {0}".format(amoptd['theseus_exe'])
        ample_exit.exit(msg)
    #
    # Scwrl
    #
    if amoptd['use_scwrl']:
        if not amoptd['scwrl_exe']:
            if sys.platform.startswith("win"):
                amoptd['scwrl_exe'] = 'Scwrl4.exe'
            else:
                amoptd['scwrl_exe'] = 'Scwrl4'
        try:
            amoptd['scwrl_exe'] = ample_util.find_exe(amoptd['scwrl_exe'])
        except Exception:
            msg = "Cannot find Scwrl executable: {0}".format(amoptd['scwrl_exe'])
            ample_exit.exit(msg)
    
    #
    # We use shelxe by default so if we can't find it we just warn and set use_shelxe to False
    #
    if amoptd['use_shelxe']:
        if not amoptd['shelxe_exe']:
            if sys.platform.startswith("win"):
                amoptd['shelxe_exe'] = 'shelxe.exe'
            else:
                amoptd['shelxe_exe'] = 'shelxe'
        try:
            amoptd['shelxe_exe'] = ample_util.find_exe(amoptd['shelxe_exe'])
        except Exception:
            msg = """*** Cannot find shelxe executable in PATH - turning off use of SHELXE. ***
    SHELXE is recommended for the best chance of success. We recommend you install shelxe from:
    http://shelx.uni-ac.gwdg.de/SHELX/
    and install it in your PATH so that AMPLE can use it.
    """
            logger.warn(msg)
            amoptd['use_shelxe'] = False
    #
    # If shelxe_rebuild is set we need use_shelxe to be set
    #
    if amoptd['shelxe_rebuild'] and not amoptd['use_shelxe']:
        msg = 'shelxe_rebuild is set but use_shelxe is False. Please make sure you have shelxe installed.'
        ample_exit.exit(msg)
    
    # Create the rosetta modeller - this runs all the checks required
    rosetta_modeller = None
    if amoptd['make_models'] or amoptd['make_frags'] or amoptd['nmr_remodel']:  # only need Rosetta if making models
        logger.info('Using ROSETTA so checking options')
        try:
            rosetta_modeller = rosetta_model.RosettaModel(optd=amoptd)
        except Exception, e:
            msg = "Error setting ROSETTA options: {0}".format(e)
            ample_exit.exit(msg)
    
    if amoptd['make_frags']:
        if amoptd['use_homs']:
            logger.info('Making fragments (including homologues)')
        else:
            logger.info('Making fragments EXCLUDING HOMOLOGUES')
    else:
        logger.info('NOT making Fragments')
    
    if amoptd['make_models']:
        logger.info('\nMaking Rosetta Models')
    else:
        logger.info('NOT making Rosetta Models')
        
        # Print out what is being done
    if amoptd['use_buccaneer']:
        logger.info('Rebuilding in Bucaneer')
    else:
        logger.info('Not rebuilding in Bucaneer')
    
    if amoptd['use_arpwarp']:
        logger.info('Rebuilding in ARP/wARP')
    else:
        logger.info('Not rebuilding in ARP/wARP')
    
    # cluster queueing
    if amoptd['submit_cluster'] and not amoptd['submit_qtype']:
        msg = 'Must use -submit_qtype argument to specify queueing system (e.g. QSUB, LSF ) if submitting to a cluster.'
        ample_exit.exit(msg)
    
    if amoptd['purge']:
        logger.info('*** Purge mode specified - all intermediate files will be deleted ***')
    
    return rosetta_modeller

def setup_ccp4(amoptd):
     # Make sure CCP4 is around
    if not "CCP4" in os.environ:
        msg = "Cannot find CCP4 installation - please make sure CCP4 is installed and the setup scripts have been run!"
        ample_exit.exit(msg)
        
    if not "CCP4_SCR" in os.environ:
        msg = "$CCP4_SCR environement variable not set - please make sure CCP4 is installed and the setup scripts have been run!"
        ample_exit.exit(msg)
        
    if not os.path.isdir(os.environ['CCP4_SCR']):
        msg = "Cannot find the $CCP4_SCR directory: {0}\nPlease make sure CCP4 is installed and the setup scripts have been run!".format(os.environ['CCP4_SCR'])
        ample_exit.exit(msg)

    # Record the CCP4 version we're running with  - also required in pyrvapi_results
    amoptd['ccp4_version'] = ample_util.ccp4_version()
    return

def main():
    """Main AMPLE routine.
    
    We require this as the multiprocessing module (only on **!!*%$$!! Windoze) requires that the main module
    can be imported. We there need ample to be a python script that can be imported, hence the main routine with
    its calling protected by the if __name__=="__main__":...
    """ 
    amopt, orig_argv = process_command_line()
    
    # Make a work directory and go there - this way all output goes into this directory
    if not os.path.exists(amopt.d['run_dir']):
        print 'Cannot find run directory: {0}'.format(amopt.d['run_dir'])
        sys.exit()
    
    print 'Making a Run Directory: checking for previous runs\n'
    amopt.d['work_dir'] = ample_util.make_workdir(amopt.d['run_dir'], ccp4_jobid=amopt.d['ccp4_jobid'])
    os.chdir(amopt.d['work_dir'])
    
    # Set up logging
    ample_log = os.path.join(amopt.d['work_dir'], 'AMPLE.log')
    amopt.d['ample_log'] = ample_log
    logger = ample_util.setup_logging(ample_log)
    
    setup_ccp4(amopt.d)
    
    # Print out Version and invocation
    logger.info(ample_util.header)
    logger.info("AMPLE version: {0}".format(version.__version__))
    logger.info("Running with CCP4 version: {0}".format(".".join([str(x) for x in amopt.d['ccp4_version']])))
    logger.info("Job started at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
    logger.info("Running on host: {0}".format(platform.node()))
    logger.info("Invoked with command-line:\n{0}\n".format(orig_argv))
    logger.info("Running in directory: {0}\n".format(amopt.d['work_dir']))
    
    # Display pyrvapi results
    pyrvapi_results.display_results(amopt.d)
    
    # Bit clunky but the rosetta_modeller object checks the rosetta options so we create it and return it if needed
    rosetta_modeller = process_options(amopt.d, logger)
    
    # Bail and clean up if we were only checking the options
    if amopt.d['dry_run']:
        logger.info('Dry run finished checking options - cleaning up...')
        os.chdir(amopt.d['run_dir'])
        shutil.rmtree(amopt.d['work_dir'])
        sys.exit(0)
    
    logger.info('All needed programs are found, continuing Run')
    
    # params used
    with open(os.path.join(amopt.d['work_dir'], 'params_used.txt'), "w") as f:
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

    # Create function for monitoring jobs - static function decorator?
    if pyrvapi_results.pyrvapi:
        def monitor():
            return pyrvapi_results.display_results(amopt.d)
    else:
        monitor = None

    # Make Rosetta fragments
    if amopt.d['make_frags']:
        rosetta_modeller.generate_fragments(amopt.d)
        amopt.d['frags_3mers'] = rosetta_modeller.frags_3mers
        amopt.d['frags_9mers'] = rosetta_modeller.frags_9mers
    
    # if NMR process models first
    # break here for NMR (frags needed but not modelling
    if amopt.d['nmr_model_in']:
        if not amopt.d['nmr_remodel']:
            if not os.path.isdir(amopt.d['models_dir']): os.mkdir(amopt.d['models_dir'])
            pdb_edit.split_pdb(amopt.d['nmr_model_in'], amopt.d['models_dir'])
            nmr.standardise_lengths(amopt.d['models_dir'])
        elif amopt.d['nmr_remodel']:
            try:
                rosetta_modeller.nmr_remodel(nmr_model_in=amopt.d['nmr_model_in'],
                                             ntimes=amopt.d['nmr_process'],
                                             alignment_file=amopt.d['alignment_file'],
                                             remodel_fasta=amopt.d['nmr_remodel_fasta'],
                                             monitor=monitor)
            except Exception, e:
                ex_type, ex, tb = sys.exc_info()
                msg = "Error remodelling NMR ensemble: {0}".format(e)
                ample_exit.exit(msg, tb)
    elif amopt.d['make_models']:
        # Make the models
        logger.info('----- making Rosetta models--------')
        logger.info('making {0} models...'.format(amopt.d['nmodels']))
        try:
            rosetta_modeller.ab_initio_model(monitor=monitor)
        except Exception, e:
            ex_type, ex, tb = sys.exc_info()
            msg = "Error running ROSETTA to create models: {0}".format(e)
            ample_exit.exit(msg, tb)
        if not pdb_edit.check_pdb_directory(amopt.d['models_dir'], sequence=amopt.d['sequence']):
            msg = "Problem with rosetta pdb files - please check the log for more information"
            ample_exit.exit(msg)
            
        msg = 'Modelling complete - models stored in: {0}\n'.format(amopt.d['models_dir'])
    elif amopt.d['import_models']:
        logger.info('Importing models from directory: {0}\n'.format(amopt.d['models_dir']))
        if amopt.d['use_scwrl']:
            msg = "Processing sidechains of imported models from {0} with Scwl\n".format(amopt.d['models_dir'])
            models_dir_scwrl = os.path.join(amopt.d['work_dir'], os.path.basename(amopt.d['models_dir']) + "_scwrl")
            if os.path.isdir(models_dir_scwrl):
                msg = "Scwrl models directory {0} already exists-please move it aside".format(models_dir_scwrl)
                ample_exit.exit(msg)
            os.mkdir(models_dir_scwrl)
            msg += "Scwrl-processed models will be placed in directory: {0}".format(models_dir_scwrl)
            msg += "Running Scwrl..."
            logger.info(msg)
            scwrl = add_sidechains_SCWRL.Scwrl(scwrlExe=amopt.d['scwrl_exe'])
            scwrl.processDirectory(inDirectory=amopt.d['models_dir'], outDirectory=models_dir_scwrl)
            amopt.d['models_dir'] = models_dir_scwrl
            logger.info("Finished processing models with Scwrl")
    
    # Do the clustering
    ensembles = []  # List of ensembles - 1 per cluster
    ensemble_options = {}
    if amopt.d['import_ensembles']:
        # Importing pre-made ensembles
        # Set list of ensembles to the one we are importing
        msg = "Importing ensembles from directory: {0}".format(amopt.d['ensembles_dir'])
        logger.info(msg)
        ensembles = glob.glob(os.path.join(amopt.d['ensembles_dir'], '*.pdb'))
        amopt.d['ensembles'] = ensembles
        amopt.d['ensembles_data'] = [ {'name' : os.path.splitext(os.path.basename(e))[0], 'ensemble_pdb' : e} for e in ensembles ]
    elif amopt.d['ideal_helices']:
        ensembles, ensemble_options, ensembles_data = ample_util.ideal_helices(amopt.d['fasta_length'])
        amopt.d['ensembles_data'] = ensembles_data
        amopt.d['ensembles'] = ensembles
        logger.info("*** Using ideal helices to solve structure ***")
    else:
        # Check we have some models to work with
        if not amopt.d['import_cluster'] and not glob.glob(os.path.join(amopt.d['models_dir'], "*.pdb")):
            ample_util.saveAmoptd(amopt.d)
            msg = "ERROR! Cannot find any pdb files in: {0}".format(amopt.d['models_dir'])
            ample_exit.exit(msg)
        
        amopt.d['ensemble_ok'] = os.path.join(amopt.d['work_dir'],'ensemble.ok')
        if amopt.d['submit_cluster']:
            # Pickle dictionary so it can be opened by the job to get the parameters
            ample_util.saveAmoptd(amopt.d)
            script = ensemble.cluster_script(amopt.d)
            ok = workers.run_scripts(job_scripts=[script],
                                     monitor=monitor,
                                     chdir=True,
                                     nproc=amopt.d['nproc'],
                                     job_time=3600,
                                     job_name='ensemble',
                                     submit_cluster=amopt.d['submit_cluster'],
                                     submit_qtype=amopt.d['submit_qtype'],
                                     submit_queue=amopt.d['submit_queue'],
                                     submit_array=amopt.d['submit_array'],
                                     submit_max_array=amopt.d['submit_max_array'])
            # queue finished so unpickle results
            with open(amopt.d['results_path'], "r") as f: amopt.d = cPickle.load(f)
        else:
            try: ensemble.create_ensembles(amopt.d)
            except Exception, e:
                msg = "Error creating ensembles: {0}".format(e)
                ample_exit.exit(msg)
                
        # Check we have something to work with
        if not os.path.isfile(amopt.d['ensemble_ok']) or not amopt.d.has_key('ensembles') or not len(amopt.d['ensembles']):
            msg = "Problem generating ensembles!"
            ample_exit.exit(msg)
            
        ensembles = amopt.d['ensembles']
        if not amopt.d['homologs']:
            ensemble_summary = ensemble.ensemble_summary(amopt.d['ensembles_data'])
            logger.info(ensemble_summary)
        
        # Save truncation levels for analysis
        
    
    # Update results view
    pyrvapi_results.display_results(amopt.d)
    
    # Save the results
    ample_util.saveAmoptd(amopt.d)
    
    #
    # Bail here if we didn't create anything
    #
    if not len(ensembles):
        msg = "### AMPLE FAILED TO GENERATE ANY ENSEMBLES! ###\nExiting..."
        ample_exit.exit(msg)
    
    # MRBUMP analysis of the ensembles
    logger.info('----- Running MRBUMP on ensembles--------\n\n')
    if len(ensembles) < 1:
        msg = "ERROR! Cannot run MRBUMP as there are no ensembles!"
        ample_exit.exit(msg)
    
    bump_dir = os.path.join(amopt.d['work_dir'], 'MRBUMP')
    if not os.path.exists(bump_dir): os.mkdir(bump_dir)
    os.chdir(bump_dir)
    amopt.d['mrbump_dir'] = bump_dir
    amopt.d['mrbump_results'] = []
    logger.info("Running MRBUMP jobs in directory: {0}".format(bump_dir))
    
    # Create job scripts
    logger.info("Generating MRBUMP runscripts")
    mrbump_jobtime = 172800  # allow 48 hours for each mrbump job
    job_scripts = mrbump_ensemble.write_mrbump_files(ensembles,
                                                     amopt.d,
                                                     job_time=mrbump_jobtime,
                                                     ensemble_options=ensemble_options)
    #ample_exit.exit("NOOOOOOO!!!!")
    
    # Create function for monitoring jobs - static function decorator?
    if pyrvapi_results.pyrvapi:
        def monitor():
            r = mrbump_results.ResultsSummary()
            r.extractResults(amopt.d['mrbump_dir'], purge=amopt.d['purge'])
            amopt.d['mrbump_results'] = r.results
            return pyrvapi_results.display_results(amopt.d)
    else:
        monitor = None
    
    ok = workers.run_scripts(job_scripts=job_scripts,
                             monitor=monitor,
                             check_success=mrbump_results.checkSuccess,
                             early_terminate=amopt.d['early_terminate'],
                             chdir=False,
                             nproc=amopt.d['nproc'],
                             job_time=mrbump_jobtime,
                             job_name='mrbump',
                             submit_cluster=amopt.d['submit_cluster'],
                             submit_qtype=amopt.d['submit_qtype'],
                             submit_queue=amopt.d['submit_queue'],
                             submit_array=amopt.d['submit_array'],
                             submit_max_array=amopt.d['submit_max_array'])

    if not ok:
        msg = "Error running MRBUMP on the ensembles!\nCheck logs in directory: {0}".format(amopt.d['mrbump_dir'])
        ample_exit.exit(msg)
        
    # Collect the MRBUMP results
    results_summary = mrbump_results.ResultsSummary()
    results_summary.extractResults(bump_dir, purge=amopt.d['purge'])
    amopt.d['mrbump_results'] = results_summary.results
    ample_util.saveAmoptd(amopt.d)
    
    # Timing data
    time_stop = time.time()
    elapsed_time = time_stop - time_start
    run_in_min = elapsed_time / 60
    run_in_hours = run_in_min / 60
    msg = '\nMR and shelx DONE\n\n ALL DONE  (in {0} hours) \n----------------------------------------\n'.format(run_in_hours)
    logging.info(msg)
    
    # Benchmark mode
    if amopt.d['benchmark_mode']:
        benchmark.analyse(amopt.d)
        ample_util.saveAmoptd(amopt.d)
        
    # Now print out the final summary
    summary = mrbump_results.finalSummary(amopt.d)
    logger.info(summary)
    
    # Finally update pyrvapi results
    pyrvapi_results.display_results(amopt.d)

    logger.info("AMPLE finished at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
    return

if __name__ == "__main__":
    main()

