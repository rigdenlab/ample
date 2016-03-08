#!/usr/bin/env ccp4-python
"""
This is AMPLE

This script is named ample.py due to a problem with running the multiprocessing (which is used to parallelise
the running of jobs on a local machine - see python/workers.py) module under windows.
The multiprocessing module on windows requires that it can import the main module, and the import machinery
requires that any file being imported is named <foo>.py, and any changes to this would require hacking the 
multiprocessing module, so to avoid this, our script must be called ample.py
"""

# python imports
import argparse
import cPickle
import glob
import logging
import os
import platform
import shutil
import sys
import time

# Our imports
import ample_config
import ample_contacts
import ample_ensemble
import ample_mrbump
import ample_exit
import ample_options
import ample_sequence
import ample_util
import ample_benchmark
import ensemble
import mtz_util
import pdb_edit
import pyrvapi_results
import rosetta_model
import version
import workers

def setup_console_logging():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # First create console logger for outputting stuff
    # create file handler and set level to debug
    # Seems they changed the api in python 2.6->2.7
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(strm=sys.stdout)
    cl.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s\n') # Always add a blank line after every print
    cl.setFormatter(formatter)
    logger.addHandler(cl)
    return logger

def setup_file_logging(main_logfile, debug_logfile):
    """
    Set up the various log files/console logging
    and return the logger
    """
    logger = logging.getLogger()

    # create file handler for debug output
    fl = logging.FileHandler(debug_logfile)
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)

    # Finally create the main logger
    fl = logging.FileHandler(main_logfile)
    fl.setLevel(logging.INFO)
    fl.setFormatter(formatter) # Same formatter as screen
    logger.addHandler(fl)
    
    return logger

logger = setup_console_logging()

class Ample(object):
    """Class to generate ensembles from ab inito models (all models must have same sequence)
    
    """
    def __init__(self):
        self.amopt = None
        self.output_gui = None
        return

    def check_mandatory_options(self, amoptd):
        """We check there here rather then with argparse as there doesn't seem to be an easy way to get the logic to work
        of having overlapping required and mutually exclusive options"""
        
        if not (amoptd['fasta'] or amoptd['restart_pkl']):
            msg = "One of -fasta  or -restart_pkl option is required."
            ample_exit.exit_error(msg)
            
        if (amoptd['contact_file'] or amoptd['bbcontacts_file']) and amoptd['restraints_file']:
            msg = "Only one option of -contact_file or -restraints_file allowed."
            ample_exit.exit_error(msg)
        
        if not amoptd['restart_pkl'] and not (amoptd['mtz'] or amoptd['sf_cif']):
            msg = "A crystallographic data file must be supplied with the -mtz or -sc_cif options."
            ample_exit.exit_error(msg)
            
        if amoptd['mtz'] and amoptd['sf_cif']:
            msg = "Please supply a single crystallographic data file."
            ample_exit.exit_error(msg)
            
        if amoptd['devel_mode'] and amoptd['quick_mode']:
            msg = "Only one of quick_mode or devel_mode is permitted"
            ample_exit.exit_error(msg)
            
        if amoptd['molrep_only'] and amoptd['phaser_only']:
            msg = "Only one of molrep_only or phaser_only is permitted"
            ample_exit.exit_error(msg)
        return
            
    def process_command_line(self, args=None):
        """Process the command-line.
        args is an optional argument that can hold the command-line arguments if we have been called from within
        python for testing.
        """
        # get command line options
        parser = argparse.ArgumentParser(description='AMPLE: Ab initio Modelling of Proteins for moLEcular replacement', prefix_chars="-")
        
        parser.add_argument('-alignment_file', type=str, nargs=1,
                           help='Alignment file in fasta format. For homologues the first line of each sequence must be the pdb file name')
        
        parser.add_argument('-allow_his_tag', metavar='True/False', type=str, nargs=1,
                           help='Allow HIS tags in the input sequence')
        
        parser.add_argument('-blast_dir', type=str, nargs=1,
                           help='Directory where ncbi blast is installed (binaries in expected in bin subdirectory)')
        
        parser.add_argument('-cluster_dir', type=str, nargs=1,
                           help='Path to directory of pre-clustered models to import')
        
        parser.add_argument('-cluster_method', type=str, nargs=1,
                           help='How to cluster the models for ensembling (spicker|fast_protein_cluster')
        
        parser.add_argument('-ccp4_jobid', type=int, nargs=1,
                           help='Set the CCP4 job id - only needed when running from the CCP4 GUI')
        
        parser.add_argument('-debug', metavar='True/False', type=str, nargs=1,
                           help='Run in debug mode (CURRENTLY UNUSED)')
    
        parser.add_argument('-devel_mode', metavar='devel_mode', type=str, nargs=1,
                           help='Preset options to run in development mode - takes longer')
        
        parser.add_argument('-domain_all_chains_pdb', type=str, nargs=1,
                           help='Fixed input to mr bump')
        
        parser.add_argument('-domain_termini_distance', type=str, nargs=1,
                           help='distance between termini for insert domains')
        
        parser.add_argument('-dry_run', metavar='True/False', type=str, nargs=1,
                             help='Check if input files and supplied options are valid.')
        
        parser.add_argument('-early_terminate', metavar='True/False', type=str, nargs=1,
                             help='Stop the run as soon as a success has been found.')
        
        parser.add_argument('-ensembles', type=str, nargs=1,
                           help='Path to directory containing existing ensembles')
        
        parser.add_argument('-fasta', type=str, nargs=1,
                           help='protein fasta file. (required)')
        
        parser.add_argument('-fast_protein_cluster_exe', type=str, nargs=1,
                           help='path to fast_protein_cluster executable')
        
        parser.add_argument('-F', metavar='flag for F', type=str, nargs=1,
                           help='Flag for F column in the MTZ file')
        
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
    
        parser.add_argument('-mr_sequence', type=str, nargs=1,
                           help="sequence file for crystal content (if different from what's given by -fasta)")
    
        parser.add_argument('-mtz', metavar='MTZ in', type=str, nargs=1,
                           help='The MTZ file with the reflection data.')
    
        parser.add_argument('-mustang_exe', metavar='mustang_exe', type=str, nargs=1,
                           help='Path to the mustang executable')
    
        parser.add_argument('-name', metavar='job_name', type=str, nargs=1,
                           help='4-letter identifier for job [ampl]')
        
        parser.add_argument('-native_pdb', metavar='native_pdb', type=str, nargs=1,
                           help='Path to the crystal structure PDB for benchmarking.')
        
        parser.add_argument('-nmodels', metavar='number of models', type=int, nargs=1,
                           help='number of models to make (default: 1000)')
        
        parser.add_argument('-nr', metavar='nr', type=str, nargs=1,
                           help='Path to the NR non-redundant sequence database')
        
        parser.add_argument('-nmr_model_in', metavar='nmr_model_in', type=str, nargs=1,
                           help='PDB with NMR models')
        
        parser.add_argument('-nmr_process', type=int, nargs=1,
                           help='number of times to process the NMR models')
        
        parser.add_argument('-nmr_remodel', metavar='True/False', type=str, nargs=1,
                           help='Remodel the NMR structures')
        
        parser.add_argument('-nmr_remodel_fasta', type=str, nargs=1,
                           help='The FASTA sequence to be used for remodelling the NMR ensemble if different from the default FASTA sequence')
        
        parser.add_argument('-no_gui', metavar='True/False', type=str, nargs=1,
                           help='Do not display the AMPLE gui.')
        
        parser.add_argument('-nproc', type=int, nargs=1,
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
    
        parser.add_argument('-psipred_ss2', metavar='PSIPRED_FILE', type=str, nargs=1,
                           help='Psipred secondary structure prediction file')
        
        parser.add_argument('-phenix_exe', type=str, nargs=1,
                           help='Path to Phenix executable')
    
        parser.add_argument('-quick_mode', metavar='True/False', type=str, nargs=1,
                           help='Preset options to run quickly, but less thoroughly')
        
        parser.add_argument('-restart_pkl', type=str, nargs=1,
                           help='Rerun a job using the pickled ample dictionary')
        
        parser.add_argument('-run_dir', metavar='run_directory', type=str, nargs=1,
                           help='Directory where the AMPLE work directory will be created [current dir]')
        
        parser.add_argument('-scwrl_exe', metavar='path to scwrl', type=str, nargs=1,
                           help='Path to Scwrl4 executable')
        
        parser.add_argument('-score_matrix', type=str, nargs=1,
                           help='Path to score matrix for spicker')
        
        parser.add_argument('-score_matrix_file_list', type=str, nargs=1,
                           help='File with list of ordered model names for the score_matrix')
    
        parser.add_argument('-sf_cif', type=str, nargs=1,
                           help='Path to a structure factor CIF file (instead of MTZ file)')
    
        parser.add_argument('-side_chain_treatments', type=str, nargs='+', action='append',
                           help='The side chain treatments to use. Default: {0}'.format(ample_ensemble.SIDE_CHAIN_TREATMENTS))
        
        parser.add_argument('-SIGF', type=str, nargs=1,
                           help='Flag for SIGF column in the MTZ file')
        
        parser.add_argument('-subcluster_program', type=str, nargs=1,
                           help='Program for subclustering models [maxcluster]')
        
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
        
        parser.add_argument('-truncation_method', type=str, nargs=1,
                           help='How to truncate the models for ensembling percent|thresh|focussed')
        
        parser.add_argument('-truncation_pruning', type=str, nargs=1,
                           help='Whether to remove isolated residues (single)')
    
        parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__))
        
        parser.add_argument('-webserver_uri', type=str, nargs=1,
                           help='URI of the webserver directory - also indicates we are running as a webserver')
        
        parser.add_argument('-work_dir', type=str, nargs=1,
                           help='Path to the directory where AMPLE will run (will be created if it doesn\'t exist)')
        
        # Contact options
        contact_group = parser.add_argument_group("Contact Restraints Options")
        
        contact_group.add_argument('-bbcontacts_file', type=str, nargs=1,
                           help='Additional bbcontacts file. Requires normal contactfile')
        
        contact_group.add_argument('-contact_file', type=str, nargs=1,
                           help='Residue contact file in CASP RR format')
        
        contact_group.add_argument('-disulfide_constraints_file', type=str, nargs=1,
                           help='Disulfide residue constraints for ab initio modelling')
        
        contact_group.add_argument('-distance_to_neighbour', type=int, nargs=1,
                           help="Min. distance between residue pairs for contact (default=5)")
        
        contact_group.add_argument('-energy_function', type=str, nargs=1,
                           help='Rosetta energy function for contact restraint conversion (default=FADE)')
        
        contact_group.add_argument('-native_cutoff', type=float, nargs=1,
                           help='Distance cutoff for reference contacts in native structure (default=8A)')
        
        contact_group.add_argument('-restraints_factor', type=float, nargs=1,
                           help='Factor (* Sequence length) determining number of contact restraints to use (default=1.0)')
        
        contact_group.add_argument('-restraints_file', type=str, nargs=1,
                           help='Residue restraints for ab initio modelling')
        
        contact_group.add_argument('-restraints_weight', type=float, nargs=1,
                           help="Additional energy weighting of restraints in Rosetta")
        
        
        # MR options
        mr_group = parser.add_argument_group('MRBUMP/Molecular Replacement Options')
        
        mr_group.add_argument('-arpwarp_cycles', type=int, nargs=1,
                           help='The number of ArpWarp cycles to run') 
        
        mr_group.add_argument('-buccaneer_cycles', type=int, nargs=1,
                           help='The number of Bucanner rebuilding cycles to run')
    
        mr_group.add_argument('-do_mr', type=str, metavar='True/False', nargs=1,
                           help='Run or skip the Molecular Replacement step')
    
        mr_group.add_argument('-molrep_only', metavar='True/False', type=str, nargs=1,
                           help='Only use Molrep for Molecular Replacement step in MRBUMP')
        
        mr_group.add_argument('-mrbump_dir', type=str, nargs=1,
                           help='Path to a directory of MRBUMP jobs (see restart_pkl)')
        
        mr_group.add_argument('-mr_keys', type=str, nargs='+', action='append',
                           help='Additional keywords for MRBUMP - are passed through without editing')
        
        mr_group.add_argument('-mr_sg_all', metavar='True/False', type=str, nargs=1,
                           help='Try all possible space groups in PHASER Molecular Replacement step in MRBUMP')
    
        mr_group.add_argument('-nmasu', type=int, nargs=1,
                           help='Manually specify the number of molecules in the asymmetric unit - sets the NMASu MRBUMP flag')
        
        mr_group.add_argument('-phaser_kill', metavar='phaser_kill', type=int, nargs=1,
                           help='Time in minutes after which phaser will be killed (0 to leave running)')
    
        mr_group.add_argument('-phaser_only', metavar='True/False', type=str, nargs=1,
                           help='Only use Phaser for Molecular Replacement step in MRBUMP')
    
        mr_group.add_argument('-phaser_rms', metavar='phaser_rms', type=float, nargs=1,
                           help='rms value for phaser (default=0.1)')
    
        mr_group.add_argument('-shelx_cycles', type=str, nargs=1,
                             help='The number of shelx cycles to run when rebuilding.')
        
        mr_group.add_argument('-shelxe_exe', metavar='path to shelxe executable', type=str, nargs=1,
                           help='Path to the shelxe executable')
        
        mr_group.add_argument('-shelxe_rebuild', metavar='True/False', type=str, nargs=1,
                           help='Rebuild shelxe traced pdb with buccaneer and arpwarp')
        
        mr_group.add_argument('-shelxe_rebuild_arpwarp', metavar='True/False', type=str, nargs=1,
                           help='Rebuild shelxe traced pdb with arpwarp')
        
        mr_group.add_argument('-shelxe_rebuild_buccaneer', metavar='True/False', type=str, nargs=1,
                           help='Rebuild shelxe traced pdb with buccaneer')
    
        mr_group.add_argument('-use_arpwarp', metavar='True/False', type=str, nargs=1,
                           help='True to use arpwarp to rebuild.')
        
        mr_group.add_argument('-use_buccaneer', metavar='True/False', type=str, nargs=1,
                           help='True to use Buccaneer')
        
        mr_group.add_argument('-use_scwrl', metavar='True/False', type=str, nargs=1,
                           help='Remodel sidechains of the decoy models using Scwrl4')
        
        mr_group.add_argument('-use_shelxe', metavar='True/False', type=str, nargs=1,
                           help='True to use shelxe')
    
        # Rosetta options
        rosetta_group = parser.add_argument_group('ROSETTA Modelling Options')
    
        rosetta_group.add_argument('-all_atom', metavar='True/False', type=str, nargs=1,
                           help="Do all-atom Rosetta modelling (adds \"-return_full_atom true\" to rosetta arguments")
        
        rosetta_group.add_argument('-frags_3mers', type=str, nargs=1,
                           help='Path to file with pre-existing Rosetta 3mer fragments')
        
        rosetta_group.add_argument('-frags_9mers', type=str, nargs=1,
                           help='Path to file with pre-existing Rosetta 3mer fragments')
        
        rosetta_group.add_argument('-make_frags', metavar='True/False', type=str, nargs=1,
                           help='set True to generate Rosetta 3mers and 9mers locally, False to import fragments')        
        
        rosetta_group.add_argument('-rg_reweight', metavar='radius of gyration reweight', type=float, nargs=1,
                           help='Set the Rosetta -rg_reweight flag to specify the radius of gyration reweight.')
        
        rosetta_group.add_argument('-rosetta_AbinitioRelax', type=str, nargs=1,
                           help='Path to Rosetta AbinitioRelax executable')
        
        rosetta_group.add_argument('-ROSETTA_cluster', type=str, nargs=1,
                           help='location of rosetta cluster')
        
        rosetta_group.add_argument('-rosetta_db', type=str, nargs=1,
                           help='Path to the Rosetta database directory')
        
        rosetta_group.add_argument('-rosetta_dir', type=str, nargs=1,
                           help='The Rosetta install directory')
    
        rosetta_group.add_argument('-rosetta_fragments_exe', type=str, nargs=1,
                           help='Location of the Rosetta make_fragments.pl script')
        
        rosetta_group.add_argument('-rosetta_version', type=float, nargs=1,
                           help='The version number of Rosetta')
    
        rosetta_group.add_argument('-transmembrane', type=str, nargs=1,
                           help='Do Rosetta modelling for transmembrane proteins')
        
        rosetta_group.add_argument('-transmembrane2', type=str, nargs=1,
                           help='Do Rosetta modelling for transmembrane proteins (NEW PROTOCOL)')
        
        rosetta_group.add_argument('-transmembrane_octopusfile', type=str, nargs=1,
                           help='Octopus transmembrane topology predicition file')
        
        rosetta_group.add_argument('-transmembrane_spanfile', type=str, nargs=1,
                           help='Span file for modelling transmembrane proteins')
        
        rosetta_group.add_argument('-transmembrane_lipofile', type=str, nargs=1,
                           help='Lips4 file for modelling transmembrane proteins')
    
        rosetta_group.add_argument('-use_homs', metavar='True/False', type=str, nargs=1,
                           help='Select ROSETTA fragments from homologous models')
    
        # convert args to dictionary
        argso = parser.parse_args(args)
        
        # Create amopt object to hold options
        amopt = ample_options.AmpleOptions()
        
        # Now put them in the amopt object - this also sets/checks any defaults
        amopt.populate(argso)
        
        return amopt
    
    def process_options(self, amoptd):
        
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
                ample_exit.exit_error(msg)
        else:
            amoptd['mr_sequence'] = amoptd['fasta']
            
        # Process the fasta file and run all the checks on the sequence    
        ample_sequence.process_fasta(amoptd)
    
        #
        # Not sure if name actually required - see make_fragments.pl
        #
        if amoptd['name'] and len(amoptd['name']) != 4:
            msg = '-name argument is the wrong length, use 4 chars eg ABCD'
            ample_exit.exit_error(msg)
            
        # Underscore required by rosetta make_fragments.pl
        amoptd['name'] += '_'
        
        ###############################################################################
        #
        # Contact file processing
        #
        ###############################################################################
        
        if amoptd['contact_file'] or amoptd['bbcontacts_file']:
            ample_contacts.checkOptions(amoptd)
            amoptd['use_contacts'] = True
        ###############################################################################
        #
        # MTZ file processing
        #
        ###############################################################################
        try:
            mtz_util.processReflectionFile(amoptd)
        except Exception, e:
            msg = "Error processing reflection file: {0}".format(e)
            ample_exit.exit_error(msg, sys.exc_info()[2])
        logger.info("Using MTZ file: {0}".format(amoptd['mtz']))
        
        ###############################################################################
        #
        # Modelling and ensemble options
        #
        ###############################################################################
        
        # Set default name for modelling directory
        amoptd['models_dir'] = os.path.join(amoptd['work_dir'], "models")
        
        # Check if importing ensembles
        if amoptd['ensembles']:
            amoptd['import_ensembles'] = True # checks are made in ensembles.import_ensembles
            amoptd['make_frags'] = False
            amoptd['make_models'] = False
        elif amoptd['cluster_dir']:
            if not os.path.isdir(amoptd['cluster_dir']):
                msg = "Import cluster cannot find directory: {0}".format(amoptd['cluster_dir'])
                ample_exit.exit_error(msg)
            if not glob.glob(os.path.join(amoptd['cluster_dir'], "*.pdb")):
                msg = "Import cluster cannot find pdbs in directory: {0}".format(amoptd['cluster_dir'])
                ample_exit.exit_error(msg)
            logger.info("Importing pre-clustered models from directory: {0}\n".format(amoptd['cluster_dir']))   
            amoptd['cluster_method'] = 'import'
            amoptd['make_frags'] = False
            amoptd['make_models'] = False
        elif amoptd['ideal_helices']:
            amoptd['make_frags'] = False
            amoptd['make_models'] = False
        elif amoptd['homologs']:
            amoptd['make_frags'] = False
            amoptd['make_models'] = False
            if not  os.path.isfile(str(amoptd['alignment_file'])):
                # We need to use gesamt or mustang to do the alignment
                if amoptd['homolog_aligner'] == 'gesamt':
                    if not ample_util.is_exe(str(amoptd['gesamt_exe'])):
                        amoptd['gesamt_exe'] = os.path.join(os.environ['CCP4'],'bin','gesamt' + ample_util.EXE_EXT)
                    if not ample_util.is_exe(str(amoptd['gesamt_exe'])):
                        msg = 'Using homologs without an alignment file and cannot find gesamt_exe: {0}'.format(amoptd['gesamt_exe'])
                        ample_exit.exit_error(msg)
                elif amoptd['homolog_aligner'] == 'mustang':
                    if not ample_util.is_exe(str(amoptd['mustang_exe'])):
                        msg = 'Using homologs without an alignment file and cannot find mustang_exe: {0}'.format(amoptd['mustang_exe'])
                        ample_exit.exit_error(msg)
                else:
                    msg = 'Unknown homolog_aligner: {0}'.format(amoptd['homolog_aligner'])
                    ample_exit.exit_error(msg)
            if not os.path.isdir(str(amoptd['models'])):
                msg = "Homologs option requires a directory of pdb models to be supplied\n" + \
                "Please supply the models with the -models flag"
                ample_exit.exit_error(msg)
            amoptd['import_models'] = True
        elif amoptd['models']:
            amoptd['import_models'] = True
            amoptd['make_frags'] = False
            amoptd['make_models'] = False
            
        # Check import flags
        if amoptd['import_ensembles'] and (amoptd['import_models']):
                msg = "Cannot import both models and ensembles/clusters!"
                ample_exit.exit_error(msg)
        
        # NMR Checks
        if amoptd['nmr_model_in']:
            msg = "Using nmr_model_in file: {0}".format(amoptd['nmr_model_in'])
            logger.info(msg)
            if not os.path.isfile(amoptd['nmr_model_in']):
                msg = "nmr_model_in flag given, but cannot find file: {0}".format(amoptd['nmr_model_in'])
                ample_exit.exit_error(msg)
            if amoptd['nmr_remodel']:
                amoptd['make_models'] = True
                if amoptd['nmr_remodel_fasta']:
                    if not os.path.isfile(amoptd['nmr_remodel_fasta']):
                        msg = "Cannot find nmr_remodel_fasta file: {0}".format(amoptd['nmr_remodel_fasta'])
                        ample_exit.exit_error(msg)
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
                        ample_exit.exit_error(msg)
                    amoptd['make_frags'] = False
    
            else:
                amoptd['make_frags'] = False
                amoptd['make_models'] = False
                msg = "Running in NMR truncate only mode"
                logger.info(msg)
    
        elif amoptd['make_models']:
            if not os.path.isdir(amoptd['models_dir']): os.mkdir(amoptd['models_dir'])
            # If the user has given both fragment files we check they are ok and unset make_frags
            if amoptd['frags_3mers'] and amoptd['frags_9mers']:
                if not os.path.isfile(amoptd['frags_3mers']) or not os.path.isfile(amoptd['frags_9mers']):
                    msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(amoptd['frags_3mers'], amoptd['frags_9mers'])
                    ample_exit.exit_error(msg)
                amoptd['make_frags'] = False
            if amoptd['make_frags'] and (amoptd['frags_3mers'] or  amoptd['frags_9mers']):
                msg = "make_frags set to true, but you have given the path to the frags_3mers or frags_9mers"
                ample_exit.exit_error(msg)
        
            if not amoptd['make_frags'] and not (amoptd['frags_3mers'] and amoptd['frags_9mers']):
                msg = """*** Missing fragment files! ***
        Please supply the paths to the fragment files using the -frags_3mers and -frags_9mers flags.
        These can be generated using the Robetta server: http://robetta.bakerlab.org
        Please see the AMPLE documentation for further information."""
                ample_exit.exit_error(msg)
        
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
                ample_exit.exit_error(msg)
    
        # MR programs
        if amoptd['molrep_only']:
                amoptd['phaser_only'] = False
                #msg = 'you say you want molrep only AND phaser only, choose one or both'
                #ample_exit.exit_error(msg)
        
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
                ample_exit.exit_error(msg)
            amoptd['benchmark_mode'] = True
            logger.info("*** AMPLE running in benchmark mode ***")
            amoptd['benchmark_dir'] = os.path.join(amoptd['work_dir'], "benchmark")
    
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
        if amoptd['cluster_method'] == 'spicker' or amoptd['cluster_method'] == 'spicker_qscore' or \
            amoptd['cluster_method'] == 'spicker_tmscore':
            if not amoptd['spicker_exe']:
                amoptd['spicker_exe'] = 'spicker'  + ample_util.EXE_EXT
            try:
                amoptd['spicker_exe'] = ample_util.find_exe(amoptd['spicker_exe'])
            except Exception:
                msg = "Cannot find spicker executable: {0}".format(amoptd['spicker_exe'])
                ample_exit.exit_error(msg)
        elif amoptd['cluster_method'] == 'fast_protein_cluster':
            if not amoptd['fast_protein_cluster_exe']: amoptd['fast_protein_cluster_exe'] = 'fast_protein_cluster'
            try:
                amoptd['fast_protein_cluster_exe'] = ample_util.find_exe(amoptd['fast_protein_cluster_exe'])
            except Exception:
                msg = "Cannot find fast_protein_cluster executable: {0}".format(amoptd['fast_protein_cluster_exe'])
                ample_exit.exit_error(msg)
        elif amoptd['cluster_method'] == 'import' or \
             amoptd['cluster_method'] == 'random':
            pass
        else:
            msg = "Unrecognised cluster_method: {0}".format(amoptd['cluster_method'])
            ample_exit.exit_error(msg)
        if not amoptd['theseus_exe']:
            amoptd['theseus_exe'] = 'theseus' + ample_util.EXE_EXT
        try:
            amoptd['theseus_exe'] = ample_util.find_exe(amoptd['theseus_exe'])
        except Exception:
            msg = "Cannot find theseus executable: {0}".format(amoptd['theseus_exe'])
            ample_exit.exit_error(msg)
        #
        # SCRWL - we always check for SCRWL as if we are processing QUARK models we want to add sidechains to them
        #
        #if amoptd['use_scwrl']:
        if not amoptd['scwrl_exe']:
            amoptd['scwrl_exe'] = 'Scwrl4' + ample_util.EXE_EXT
        try:
            amoptd['scwrl_exe'] = ample_util.find_exe(amoptd['scwrl_exe'])
        except Exception as e:
            logger.info("Cannot find Scwrl executable: {0}".format(amoptd['scwrl_exe']))
            if amoptd['use_scwrl']: raise(e)
        #
        # We use shelxe by default so if we can't find it we just warn and set use_shelxe to False
        #
        if amoptd['use_shelxe']:
            if not amoptd['shelxe_exe']:
                amoptd['shelxe_exe'] = 'shelxe' + ample_util.EXE_EXT
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
            ample_exit.exit_error(msg)
        
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
            ample_exit.exit_error(msg)
        
        if amoptd['purge']:
            logger.info('*** Purge mode specified - all intermediate files will be deleted ***')
        
        return
    
    def process_restart_options(self, amoptd):
        """
        For any new command-line options, we update the old dictionary with the new values
        We then go through the new dictionary and set ant of the flags corresponding to the data we find:
        
        
        if restart.pkl
        - if completed mrbump jobs
            make_frags, make_models, make_ensembles = False
            make_mr = True
          - if all jobs aren't completed, rerun the remaining mrbump jobs - IN THE OLD DIRECTORY?
          - if all jobs are completed and we are in benchmark mode run the benchmarking
            make_frags, make_models, make_ensembles, make_mr = False
            make_benchmark = True
          - END
        - if ensemble files
           - if no ensemble data, create ensemble data
           make_frags, make_models, make_ensembles = False
           make_mr = True
           - create and run the mrbump jobs - see above
           
           # BElow all same as default
        - if models and no ensembles
          - create ensembles from the models
        
        FLAGS
        make_frags
        make_models
        make_ensembles
        make_mr
        make_benchmark
        
        We return the dictionary as we may need to change it and it seems we can't change the extermal
        reference in this scope. I think?...
        """
        if not amoptd['restart_pkl']: return amoptd
        if not os.path.isfile(amoptd['restart_pkl']):
            msg = 'Cannot find restart_pkl file: {0}'.format(amoptd['restart_pkl'])
            ample_exit.exit_error(msg)
        
        logger.info('Restarting from existing pkl file: {0}'.format(amoptd['restart_pkl']))
        # We use the old dictionary, but udpate it with any new values
        with open(amoptd['restart_pkl']) as f: amoptd_old = cPickle.load(f)
        
        # Update key variables that differ with a new run - everything else uses the old values
        amoptd_old['ample_log'] = amoptd['ample_log']
        amoptd_old['run_dir'] = amoptd['run_dir']
        amoptd_old['work_dir'] = amoptd['work_dir']
        amoptd_old['benchmark_mode'] = amoptd['benchmark_mode']
        amoptd_old['benchmark_dir'] = os.path.join(amoptd['work_dir'], "benchmark")
        amoptd_old['results_path'] = os.path.join(amoptd['work_dir'],'resultsd.pkl')
        
        # Now update any variables that were given on the command-line
        for k in amoptd['cmdline_flags']:
            logger.debug("Restart updating amopt variable: {0} : {1}".format(k, amoptd[k]))
            amoptd_old[k] = amoptd[k]
        
        # We can now replace the old dictionary with this new one
        amoptd = amoptd_old
        
        # Go through and see what we need to do
        
        # Reset all variables for doing stuff - otherwise we will always restart from the earliest point
        amoptd['make_ensembles'] = False
        amoptd['import_ensembles'] = False # Needs thinking about - have to set so we don't just reimport models/ensembles
        amoptd['import_models'] = False # Needs thinking about
        amoptd['make_models'] = False
        amoptd['make_frags'] = False
        
        # First see if we should benchmark this job. The user may not have supplied a native_pdb with the original
        # job and we only set benchmark mode on seeing the native_pdb
        if amoptd['native_pdb']:
            if not os.path.isfile(amoptd['native_pdb']):
                msg = "Cannot find native_pdb: {0}".format(amoptd['native_pdb'])
                logger.critical(msg)
                raise RuntimeError(msg)
            amoptd['benchmark_mode'] = True
            logger.info('Restart using benchmark mode')
            
        # We always check first to see if there are any mrbump jobs
        amoptd['mrbump_scripts'] = []
        if 'mrbump_dir' in amoptd:
            amoptd['mrbump_scripts'] = ample_mrbump.unfinished_scripts(amoptd)
            if not amoptd['mrbump_scripts']:
                amoptd['do_mr'] = False
    
        if amoptd['do_mr']:
            if len(amoptd['mrbump_scripts']):
                logger.info('Restarting from unfinished mrbump scripts: {0}'.format(amoptd['mrbump_scripts']))
            elif 'ensembles' in amoptd and amoptd['ensembles'] and len(amoptd['ensembles']):
                # Rerun from ensembles - check for data/ensembles are ok?
                logger.info('Restarting from existing ensembles: {0}'.format(amoptd['ensembles']))
            elif amoptd['models_dir'] and amoptd['models_dir'] and os.path.isdir(amoptd['models_dir']):
                logger.info('Restarting from existing models: {0}'.format(amoptd['models_dir']))
                # Check the models
                allsame = False if amoptd['homologs'] else True 
                if not pdb_edit.check_pdb_directory(amoptd['models_dir'], sequence=None, single=True, allsame=allsame):
                    msg = "Error importing restart models: {0}".format(amoptd['models_dir'])
                    ample_exit.exit_error(msg)
                amoptd['make_ensembles'] = True
            elif amoptd['frags_3mers'] and amoptd['frags_9mers']:
                logger.info('Restarting from existing fragments: {0}, {1}'.format(amoptd['frags_3mers'], amoptd['frags_9mers']))
                amoptd['make_models'] = True
        
        return amoptd
    
    def process_rosetta_options(self, amoptd):
        # Create the rosetta modeller - this runs all the checks required
        rosetta_modeller = None
        if amoptd['make_models'] or amoptd['make_frags']:  # only need Rosetta if making models
            logger.info('Using ROSETTA so checking options')
            try:
                rosetta_modeller = rosetta_model.RosettaModel(optd=amoptd)
            except Exception, e:
                msg = "Error setting ROSETTA options: {0}".format(e)
                ample_exit.exit_error(msg)
        return rosetta_modeller
    
    def setup_ccp4(self, amoptd):
        """Check CCP4 is available and return the top CCP4 directory"""
        # Make sure CCP4 is around
        if not "CCP4" in os.environ:
            msg = "Cannot find CCP4 installation - please make sure CCP4 is installed and the setup scripts have been run!"
            ample_exit.exit_error(msg)
            
        if not "CCP4_SCR" in os.environ:
            msg = "$CCP4_SCR environement variable not set - please make sure CCP4 is installed and the setup scripts have been run!"
            ample_exit.exit_error(msg)
            
        if not os.path.isdir(os.environ['CCP4_SCR']):
            msg = "*** WARNING ***\n"
            msg += "Cannot find the $CCP4_SCR directory: {0}\n".format(os.environ['CCP4_SCR'])
            msg += "The directory will be created, but it should have already been created by the CCP4 startup scripts\n"
            msg += "Please make sure CCP4 is installed and the setup scripts have been run."
            logger.critical(msg)
            os.mkdir(os.environ['CCP4_SCR'])
            #ample_exit.exit_error(msg)
    
        # Record the CCP4 version we're running with  - also required in pyrvapi_results
        amoptd['ccp4_version'] = ample_util.ccp4_version()
        
        return os.environ['CCP4']

    def run(self, args=None):
        """Main AMPLE routine.
        
        We require this as the multiprocessing module (only on **!!*%$$!! Windoze) requires that the main module
        can be imported. We there need ample to be a python script that can be imported, hence the main routine with
        its calling protected by the if __name__=="__main__":...
        
        args is an option argument that can contain the command-line arguments for the program - required for testing.
        """ 
        amopt = self.process_command_line(args=args)
        self.amopt = amopt
        
        # Make a work directory - this way all output goes into this directory
        if amopt.d['work_dir']:
            logger.info('Making a named work directory: {0}'.format(amopt.d['work_dir']))
            try:
                os.mkdir(amopt.d['work_dir'])
            except:
                msg = "Cannot create work_dir {0}".format(amopt.d['work_dir'])
                ample_exit.exit_error(msg, sys.exc_info()[2])
        else:
            if not os.path.exists(amopt.d['run_dir']):
                msg = 'Cannot find run directory: {0}'.format(amopt.d['run_dir'])
                ample_exit.exit_error(msg, sys.exc_info()[2])
            logger.info('Making a run directory: checking for previous runs...')
            amopt.d['work_dir'] = ample_util.make_workdir(amopt.d['run_dir'], ccp4_jobid=amopt.d['ccp4_jobid'])
        # Go to the work directory
        os.chdir(amopt.d['work_dir'])
        
        # Set up logging
        ample_log = os.path.join(amopt.d['work_dir'], 'AMPLE.log')
        debug_log = os.path.join(amopt.d['work_dir'], 'debug.log')
        amopt.d['ample_log'] = ample_log
        
        setup_file_logging(ample_log, debug_log)
        
        # Make sure the CCP4 environment is set up properly
        ccp4_home = self.setup_ccp4(amopt.d)
        
        # Print out Version and invocation
        logger.info(ample_util.header)
        logger.info("AMPLE version: {0}".format(version.__version__))
        logger.info("Running with CCP4 version: {0} from directory: {1}".format(".".join([str(x) for x in amopt.d['ccp4_version']]), ccp4_home))
        logger.info("Job started at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        logger.info("Running on host: {0}".format(platform.node()))
        logger.info("Invoked with command-line:\n{0}\n".format(" ".join(sys.argv)))
        logger.info("Running in directory: {0}\n".format(amopt.d['work_dir']))
        
        # Display pyrvapi results
        if pyrvapi_results.pyrvapi:
            self.output_gui = pyrvapi_results.AmpleOutput()
            self.output_gui.display_results(amopt.d)
        
        # Check mandatory/exclusive options
        self.check_mandatory_options(amopt.d)
        
        # Check if we are restarting from an existing pkl file - we don't process the options from this
        # run if so
        amopt.d = self.process_restart_options(amopt.d)
        if not amopt.d['restart_pkl']:
            self.process_options(amopt.d) # Only process the remaining options if we aren't in restart mode
        rosetta_modeller = self.process_rosetta_options(amopt.d)
        
        # Bail and clean up if we were only checking the options
        if amopt.d['dry_run']:
            logger.info('Dry run finished checking options - cleaning up...')
            os.chdir(amopt.d['run_dir'])
            shutil.rmtree(amopt.d['work_dir'])
            sys.exit(0)
        
        logger.info('All needed programs are found, continuing...')
        
        # Display the parameters used
        logger.debug(amopt.prettify_parameters())
        
        #######################################################
        #
        # SCRIPT PROPER STARTS HERE
        #
        ######################################################
        time_start = time.time()
    
        # Create function for monitoring jobs - static function decorator?
        if self.output_gui:
            def monitor():
                return self.output_gui.display_results(amopt.d)
        else:
            monitor = None
            
        if amopt.d['benchmark_mode']:
            # Process the native before we do anything else
            ample_benchmark.analysePdb(amopt.d)       
    
        # Make Rosetta fragments
        if amopt.d['make_frags']:
            rosetta_modeller.generate_fragments(amopt.d)
            amopt.d['frags_3mers'] = rosetta_modeller.frags_3mers
            amopt.d['frags_9mers'] = rosetta_modeller.frags_9mers
            amopt.d['psipred_ss2'] = rosetta_modeller.psipred_ss2
    
        # In case file created above we need to tell the rosetta_modeller where it is
        # otherwise not used as not created before object initialised    
        if amopt.d['make_models'] and (amopt.d['use_contacts'] or amopt.d['restraints_file']):
            cm = ample_contacts.Contacter(optd=amopt.d)
            
            cm.process_restraintsfile() if not amopt.d['use_contacts'] and amopt.d['restraints_file'] \
                else cm.process_contactfile()
    
            amopt.d['restraints_file'] = cm.restraints_file
            amopt.d['contact_map'] = cm.contact_map
            amopt.d['contact_ppv'] = cm.contact_ppv
                
        
        if amopt.d['make_models'] and amopt.d['restraints_file']: 
            rosetta_modeller.restraints_file = amopt.d['restraints_file']
        
        # if NMR process models first
        # break here for NMR (frags needed but not modelling
        if amopt.d['nmr_model_in'] and not amopt.d['nmr_remodel']:
            pdb_edit.prepare_nmr_model(amopt.d['nmr_model_in'], amopt.d['models_dir'])
        elif amopt.d['make_models']:
            # Make the models
            logger.info('----- making Rosetta models--------')
            if amopt.d['nmr_remodel']:
                try:
                    rosetta_modeller.nmr_remodel(nmr_model_in=amopt.d['nmr_model_in'],
                                                 ntimes=amopt.d['nmr_process'],
                                                 alignment_file=amopt.d['alignment_file'],
                                                 remodel_fasta=amopt.d['nmr_remodel_fasta'],
                                                 monitor=monitor)
                except Exception, e:
                    msg = "Error remodelling NMR ensemble: {0}".format(e)
                    ample_exit.exit_error(msg, sys.exc_info()[2])
            else:
                logger.info('making {0} models...'.format(amopt.d['nmodels']))
                try:
                    rosetta_modeller.ab_initio_model(monitor=monitor)
                except Exception, e:
                    msg = "Error running ROSETTA to create models: {0}".format(e)
                    ample_exit.exit_error(msg, sys.exc_info()[2])
                if not pdb_edit.check_pdb_directory(amopt.d['models_dir'], sequence=amopt.d['sequence']):
                    msg = "Problem with rosetta pdb files - please check the log for more information"
                    ample_exit.exit_error(msg)
                msg = 'Modelling complete - models stored in: {0}\n'.format(amopt.d['models_dir'])
            
        elif amopt.d['import_models']:
            logger.info('Importing models from directory: {0}\n'.format(amopt.d['models_dir']))
            if amopt.d['homologs']:
                amopt.d['models_dir'] = ample_util.extract_models(amopt.d, sequence=None, single=True, allsame=False)
            else:
                amopt.d['models_dir'] = ample_util.extract_models(amopt.d)
                # Need to check if Quark and handle things accordingly
                if amopt.d['quark_models']:
                    # We always add sidechains to QUARK models if SCWRL is installed
                    if ample_util.is_exe(amopt.d['scwrl_exe']):
                        amopt.d['use_scwrl'] = True
                    else:
                        # No SCWRL so don't do owt with the side chains
                        logger.info('Using QUARK models but SCWRL is not installed so only using {0} sidechains'.format(ample_ensemble.UNMODIFIED))
                        amopt.d['side_chain_treatments'] = [ ample_ensemble.UNMODIFIED ]
    
        # Save the results
        ample_util.saveAmoptd(amopt.d)
        
        if amopt.d['make_ensembles']:
            if amopt.d['import_ensembles']:
                ensemble.import_ensembles(amopt.d)
            elif amopt.d['ideal_helices']:
                amopt.d['ensembles'], amopt.d['ensemble_options'], amopt.d['ensembles_data'] = ample_util.ideal_helices(amopt.d['fasta_length'])
                logger.info("*** Using ideal helices to solve structure ***")
            else:
                # Check we have some models to work with
                if not amopt.d['cluster_method'] is 'import' and not glob.glob(os.path.join(amopt.d['models_dir'], "*.pdb")):
                    ample_util.saveAmoptd(amopt.d)
                    msg = "ERROR! Cannot find any pdb files in: {0}".format(amopt.d['models_dir'])
                    ample_exit.exit_error(msg)
                amopt.d['ensemble_ok'] = os.path.join(amopt.d['work_dir'],'ensemble.ok')
                if amopt.d['submit_cluster']:
                    # Pickle dictionary so it can be opened by the job to get the parameters
                    ample_util.saveAmoptd(amopt.d)
                    script = ensemble.cluster_script(amopt.d)
                    workers.run_scripts(job_scripts=[script],
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
                        ample_exit.exit_error(msg, sys.exc_info()[2])
                        
                # Check we have something to work with
                if not os.path.isfile(amopt.d['ensemble_ok']) or not amopt.d.has_key('ensembles') or not len(amopt.d['ensembles']):
                    msg = "Problem generating ensembles!"
                    ample_exit.exit_error(msg)
                    
                if not amopt.d['homologs']:
                    ensemble_summary = ensemble.ensemble_summary(amopt.d['ensembles_data'])
                    logger.info(ensemble_summary)
                
            # Save the results
            ample_util.saveAmoptd(amopt.d)
            
            # Bail here if we didn't create anything
            if not len(amopt.d['ensembles']):
                msg = "### AMPLE FAILED TO GENERATE ANY ENSEMBLES! ###\nExiting..."
                ample_exit.exit_error(msg)
        
        # Update results view
        if self.output_gui: self.output_gui.display_results(amopt.d)
         
        if amopt.d['do_mr']:
            if not amopt.d['mrbump_scripts']:
                # MRBUMP analysis of the ensembles
                logger.info('----- Running MRBUMP on ensembles--------\n\n')
                if len(amopt.d['ensembles']) < 1:
                    msg = "ERROR! Cannot run MRBUMP as there are no ensembles!"
                    ample_exit.exit_error(msg)
                 
                if amopt.d['mrbump_dir'] is None:
                    bump_dir = os.path.join(amopt.d['work_dir'], 'MRBUMP')
                    amopt.d['mrbump_dir'] = bump_dir
                else:
                    bump_dir = amopt.d['mrbump_dir']
                if not os.path.exists(bump_dir): os.mkdir(bump_dir)
                 
                amopt.d['mrbump_results'] = []
                logger.info("Running MRBUMP jobs in directory: {0}".format(bump_dir))
                
                # Sort the ensembles in a favourable way
                logger.info("Sorting ensembles")
                ensemble_pdbs_sorted = ensemble.sort_ensembles(amopt.d['ensembles'],
                                                               amopt.d['ensembles_data'])

                # Create job scripts
                logger.info("Generating MRBUMP runscripts")
                amopt.d['mrbump_scripts'] = ample_mrbump.write_mrbump_files(ensemble_pdbs_sorted,
                                                                            amopt.d,
                                                                            job_time=ample_mrbump.MRBUMP_RUNTIME,
                                                                            ensemble_options=amopt.d['ensemble_options'],
                                                                            directory=bump_dir )
            # Create function for monitoring jobs - static function decorator?
            if self.output_gui:
                def monitor():
                    r = ample_mrbump.ResultsSummary()
                    r.extractResults(amopt.d['mrbump_dir'], purge=amopt.d['purge'])
                    amopt.d['mrbump_results'] = r.results
                    return self.output_gui.display_results(amopt.d)
            else:
                monitor = None
                
            # Save results here so that we have the list of scripts and mrbump directory set
            ample_util.saveAmoptd(amopt.d)
                
            # Change to mrbump directory before running
            os.chdir(amopt.d['mrbump_dir'])  
            ok = workers.run_scripts(job_scripts=amopt.d['mrbump_scripts'],
                                     monitor=monitor,
                                     check_success=ample_mrbump.checkSuccess,
                                     early_terminate=amopt.d['early_terminate'],
                                     chdir=False,
                                     nproc=amopt.d['nproc'],
                                     job_time=ample_mrbump.MRBUMP_RUNTIME,
                                     job_name='mrbump',
                                     submit_cluster=amopt.d['submit_cluster'],
                                     submit_qtype=amopt.d['submit_qtype'],
                                     submit_queue=amopt.d['submit_queue'],
                                     submit_array=amopt.d['submit_array'],
                                     submit_max_array=amopt.d['submit_max_array'])
         
            if not ok:
                msg = "Error running MRBUMP on the ensembles!\nCheck logs in directory: {0}".format(amopt.d['mrbump_dir'])
                ample_exit.exit_error(msg)
        
            # Collect the MRBUMP results
            results_summary = ample_mrbump.ResultsSummary()
            amopt.d['mrbump_results'] = results_summary.extractResults(amopt.d['mrbump_dir'], purge=amopt.d['purge'])
            amopt.d['success'] = results_summary.success
            
            ample_util.saveAmoptd(amopt.d)
        
            # Now print out the final summary
            summary = ample_mrbump.finalSummary(amopt.d)
            logger.info(summary)
            
        # Timing data
        time_stop = time.time()
        elapsed_time = time_stop - time_start
        run_in_min = elapsed_time / 60
        run_in_hours = run_in_min / 60
        msg = os.linesep + 'All processing completed  (in {0:6.2F} hours)'.format(run_in_hours) + os.linesep
        msg += '----------------------------------------' + os.linesep
        logging.info(msg)
        
        # Benchmark mode
        if amopt.d['benchmark_mode']:
            if amopt.d['submit_cluster']:
                # Pickle dictionary so it can be opened by the job to get the parameters
                ample_util.saveAmoptd(amopt.d)
                script = ample_benchmark.cluster_script(amopt.d)
                workers.run_scripts(job_scripts=[script],
                                    monitor=monitor,
                                    chdir=True,
                                    nproc=amopt.d['nproc'],
                                    job_time=7200,
                                    job_name='benchmark',
                                    submit_cluster=amopt.d['submit_cluster'],
                                    submit_qtype=amopt.d['submit_qtype'],
                                    submit_queue=amopt.d['submit_queue'],
                                    submit_array=amopt.d['submit_array'],
                                    submit_max_array=amopt.d['submit_max_array'])
                # queue finished so unpickle results
                with open(amopt.d['results_path'], "r") as f: amopt.d = cPickle.load(f)
            else:
                ample_benchmark.analyse(amopt.d)
                ample_util.saveAmoptd(amopt.d)
        
        # Flag to show that we reached the end without error - useful for integration testing
        amopt.d['AMPLE_finished'] = True
        ample_util.saveAmoptd(amopt.d)
        
        logger.info("AMPLE finished at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        logger.info(ample_util.footer)
        
        # Finally update pyrvapi results
        if self.output_gui: self.output_gui.display_results(amopt.d)
        return

if __name__ == "__main__":
    try:
        Ample().run()
    except Exception as e:
        msg = "Error running main AMPLE program: {0}".format(e.message)
        ample_exit.exit_error(msg, sys.exc_info()[2])
