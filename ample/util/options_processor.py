"""
@author: jmht, hlfsimko
"""

import cPickle
import glob
import logging
import os
import shutil
import sys

from ample.modelling import rosetta_model
from ample.util import ample_util
from ample.util import contacts_util
from ample.util import exit_util
from ample.util import mrbump_util
from ample.util import mtz_util
from ample.util import pdb_edit
from ample.util import sequence_util

LOGGER = logging.getLogger(__name__)

def check_mandatory_options(optd):
    """We check there here rather then with argparse as there doesn't seem 
    to be an easy way to get the logic to work of having overlapping 
    required and mutually exclusive options
    """
    def _exit(msg, wdir):
        exit_util.exit_error(msg)
    
    if not (optd['fasta'] or optd['restart_pkl']):
        msg = "One of -fasta  or -restart_pkl option is required."
        _exit(msg, optd['work_dir'])
        
    if (optd['contact_file'] or optd['bbcontacts_file']) and optd['restraints_file']:
        msg = "Only one option of -contact_file or -restraints_file allowed."
        _exit(msg, optd['work_dir'])
    
    if not optd['restart_pkl'] and not (optd['mtz'] or optd['sf_cif']):
        msg = "A crystallographic data file must be supplied with the -mtz or -sc_cif options."
        _exit(msg, optd['work_dir'])
        
    if optd['do_mr'] and (optd['mtz'] and optd['sf_cif']):
        msg = "Please supply a single crystallographic data file."
        _exit(msg, optd['work_dir'])
    
    if optd['devel_mode'] and optd['quick_mode']:
        msg = "Only one of quick_mode or devel_mode is permitted"
        _exit(msg, optd['work_dir'])
        
    if optd['molrep_only'] and optd['phaser_only']:
        msg = "Only one of molrep_only or phaser_only is permitted"
        _exit(msg, optd['work_dir'])

    if optd['single_model'] and not (optd['truncation_scorefile'] and optd['truncation_scorefile_header']):
        msg = "Truncating a single model requires -truncation_scorefile and -truncation_scorefile_header"
        _exit(msg, optd['work_dir'])

    return

def process_options(optd):
    
    # Path for pickling results
    optd['results_path'] = os.path.join(optd['work_dir'], "resultsd.pkl")
    
    ###############################################################################
    #
    # FASTA processing
    #
    ###############################################################################
    # Check to see if mr_sequence was given and if not mr_sequence defaults to fasta
    if optd['mr_sequence'] != None:
        if not (os.path.exists(str(optd['mr_sequence']))):
            msg = 'Cannot find mr sequence file: {0}'.format(optd['mr_sequence'])
            exit_util.exit_error(msg)
    else:
        optd['mr_sequence'] = optd['fasta']
        
    # Process the fasta file and run all the checks on the sequence    
    sequence_util.process_fasta(optd)

    #
    # Not sure if name actually required - see make_fragments.pl
    #
    if optd['name'] and len(optd['name']) != 4:
        msg = '-name argument is the wrong length, use 4 chars eg ABCD'
        exit_util.exit_error(msg)
        
    # Underscore required by rosetta make_fragments.pl
    optd['name'] += '_'
    
    ###############################################################################
    #
    # Contact file processing
    #
    ###############################################################################
    
    if optd['contact_file'] or optd['bbcontacts_file']:
        contacts_util.checkOptions(optd)
        optd['use_contacts'] = True
    ###############################################################################
    #
    # MTZ file processing
    #
    ###############################################################################
    try:
        mtz_util.processReflectionFile(optd)
    except Exception, e:
        msg = "Error processing reflection file: {0}".format(e)
        exit_util.exit_error(msg, sys.exc_info()[2])
    
    ###############################################################################
    #
    # Modelling and ensemble options
    #
    ###############################################################################
    
    # Set default name for modelling directory
    optd['models_dir'] = os.path.join(optd['work_dir'], "models")
    
    # Check if importing ensembles
    if optd['ensembles']:
        optd['import_ensembles'] = True # checks are made in ensembles.import_ensembles
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['cluster_dir']:
        if not os.path.isdir(optd['cluster_dir']):
            msg = "Import cluster cannot find directory: {0}".format(optd['cluster_dir'])
            exit_util.exit_error(msg)
        if not glob.glob(os.path.join(optd['cluster_dir'], "*.pdb")):
            msg = "Import cluster cannot find pdbs in directory: {0}".format(optd['cluster_dir'])
            exit_util.exit_error(msg)
        LOGGER.info("Importing pre-clustered models from directory: {0}\n".format(optd['cluster_dir']))   
        optd['cluster_method'] = 'import'
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['ideal_helices']:
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['homologs']:
        optd['make_frags'] = False
        optd['make_models'] = False
        if not  os.path.isfile(str(optd['alignment_file'])):
            # We need to use gesamt or mustang to do the alignment
            if optd['homolog_aligner'] == 'gesamt':
                if not ample_util.is_exe(str(optd['gesamt_exe'])):
                    optd['gesamt_exe'] = os.path.join(os.environ['CCP4'],'bin','gesamt' + ample_util.EXE_EXT)
                if not ample_util.is_exe(str(optd['gesamt_exe'])):
                    msg = 'Using homologs without an alignment file and cannot find gesamt_exe: {0}'.format(optd['gesamt_exe'])
                    exit_util.exit_error(msg)
            elif optd['homolog_aligner'] == 'mustang':
                if not ample_util.is_exe(str(optd['mustang_exe'])):
                    msg = 'Using homologs without an alignment file and cannot find mustang_exe: {0}'.format(optd['mustang_exe'])
                    exit_util.exit_error(msg)
            else:
                msg = 'Unknown homolog_aligner: {0}'.format(optd['homolog_aligner'])
                exit_util.exit_error(msg)
        if not os.path.isdir(str(optd['models'])):
            msg = "Homologs option requires a directory of pdb models to be supplied\n" + \
            "Please supply the models with the -models flag"
            exit_util.exit_error(msg)
        optd['import_models'] = True
    elif optd['models']:
        optd['import_models'] = True
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['single_model']:
        optd['cluster_method'] = "skip"
        optd['make_frags'] = False
        optd['make_models'] = False
        optd['single_model_mode'] = True
        if optd['truncation_scorefile'] and optd['truncation_scorefile_header']:
            optd['truncation_method'] = "scores"
        
    # Check import flags
    if optd['import_ensembles'] and (optd['import_models']):
            msg = "Cannot import both models and ensembles/clusters!"
            exit_util.exit_error(msg)
    
    # NMR Checks
    if optd['nmr_model_in']:
        msg = "Using nmr_model_in file: {0}".format(optd['nmr_model_in'])
        LOGGER.info(msg)
        if not os.path.isfile(optd['nmr_model_in']):
            msg = "nmr_model_in flag given, but cannot find file: {0}".format(optd['nmr_model_in'])
            exit_util.exit_error(msg)
        if optd['nmr_remodel']:
            optd['make_models'] = True
            if optd['nmr_remodel_fasta']:
                if not os.path.isfile(optd['nmr_remodel_fasta']):
                    msg = "Cannot find nmr_remodel_fasta file: {0}".format(optd['nmr_remodel_fasta'])
                    exit_util.exit_error(msg)
            else:
                optd['nmr_remodel_fasta'] = optd['fasta']
            msg = "NMR model will be remodelled with ROSETTA using the sequence from: {0}".format(optd['nmr_remodel_fasta'])
            LOGGER.info(msg)
            
            if not optd['frags_3mers'] and optd['frags_9mers']:
                optd['make_frags'] = True
                msg = "nmr_remodel - will be making our own fragment files"
                LOGGER.info(msg)
            else:
                if not os.path.isfile(optd['frags_3mers']) or not os.path.isfile(optd['frags_9mers']):
                    msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(optd['frags_3mers'], optd['frags_9mers'])
                    exit_util.exit_error(msg)
                optd['make_frags'] = False

        else:
            optd['make_frags'] = False
            optd['make_models'] = False
            msg = "Running in NMR truncate only mode"
            LOGGER.info(msg)

    elif optd['make_models']:
        if not os.path.isdir(optd['models_dir']): os.mkdir(optd['models_dir'])
        # If the user has given both fragment files we check they are ok and unset make_frags
        if optd['frags_3mers'] and optd['frags_9mers']:
            if not os.path.isfile(optd['frags_3mers']) or not os.path.isfile(optd['frags_9mers']):
                msg = "frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(optd['frags_3mers'], optd['frags_9mers'])
                exit_util.exit_error(msg)
            optd['make_frags'] = False
        if optd['make_frags'] and (optd['frags_3mers'] or  optd['frags_9mers']):
            msg = "make_frags set to true, but you have given the path to the frags_3mers or frags_9mers"
            exit_util.exit_error(msg)
    
        if not optd['make_frags'] and not (optd['frags_3mers'] and optd['frags_9mers']):
            msg = """*** Missing fragment files! ***
    Please supply the paths to the fragment files using the -frags_3mers and -frags_9mers flags.
    These can be generated using the Robetta server: http://robetta.bakerlab.org
    Please see the AMPLE documentation for further information."""
            exit_util.exit_error(msg)
    
    ###############################################################################
    #
    # Misc options
    #
    ###############################################################################
    
    # Missing domains
    if optd['missing_domain']:
        LOGGER.info('Processing missing domain\n')
        if not os.path.exists(optd['domain_all_chains_pdb']):
            msg = 'Cannot find file domain_all_chains_pdb: {0}'.format(optd['domain_all_chains_pdb'])
            exit_util.exit_error(msg)

    # Molecular Replacement Options
    if optd['molrep_only']:
            optd['phaser_only'] = False
            #msg = 'you say you want molrep only AND phaser only, choose one or both'
            #exit_util.exit_error(msg)
    
    if optd['molrep_only']:
        optd['mrbump_programs'] = [ 'molrep' ]
    elif optd['phaser_only']:
        optd['mrbump_programs'] = [ 'phaser' ]
    else:
        optd['mrbump_programs'] = ['molrep', 'phaser']
        
    if optd['phaser_rms'] != 'auto':
        try:
            phaser_rms = float(optd['phaser_rms'])
            optd['phaser_rms'] = phaser_rms
        except ValueError as e:
            msg = "Error converting phaser_rms '{0}' to floating point: {1}".format(optd['phaser_rms'], e)
            exit_util.exit_error(msg)
    
    #
    # Benchmark Mode
    #
    if optd['native_pdb'] or optd['benchmark_mode']:
        if optd['native_pdb'] and not os.path.isfile(optd['native_pdb']):
            msg = "Cannot find crystal structure PDB: {0}".format(optd['native_pdb'])
            exit_util.exit_error(msg)
        optd['benchmark_mode'] = True
        optd['benchmark_dir'] = os.path.join(optd['work_dir'], "benchmark")
        LOGGER.info("*** AMPLE running in benchmark mode ***")
        

    ###############################################################################
    #
    # Program defaults
    #
    #
    ###############################################################################
    
    # Model building programs
    if optd['use_arpwarp']:
        if not (os.environ.has_key('warpbin') and os.path.isfile(os.path.join(os.environ['warpbin'], "auto_tracing.sh"))):
            LOGGER.warn('Cannot find arpwarp script! Disabling use of arpwarp.')
            optd['use_arpwarp'] = False
        else:
            LOGGER.info('Using arpwarp script: {0}'.format(os.path.join(os.environ['warpbin'], "auto_tracing.sh")))
    #
    # Check we can find all the required programs
    #
    # Maxcluster handled differently as we may need to download the binary
    optd['maxcluster_exe'] = ample_util.find_maxcluster(optd)
    
    #
    # Ensemble options
    #
    if optd['cluster_method'] in ['spicker', 'spicker_qscore', 'spicker_tmscore']:
        if not optd['spicker_exe']:
            optd['spicker_exe'] = 'spicker'  + ample_util.EXE_EXT
        try:
            optd['spicker_exe'] = ample_util.find_exe(optd['spicker_exe'])
        except Exception:
            msg = "Cannot find spicker executable: {0}".format(optd['spicker_exe'])
            exit_util.exit_error(msg)
    elif optd['cluster_method'] in ['fast_protein_cluster']:
        if not optd['fast_protein_cluster_exe']: optd['fast_protein_cluster_exe'] = 'fast_protein_cluster'
        try:
            optd['fast_protein_cluster_exe'] = ample_util.find_exe(optd['fast_protein_cluster_exe'])
        except Exception:
            msg = "Cannot find fast_protein_cluster executable: {0}".format(optd['fast_protein_cluster_exe'])
            exit_util.exit_error(msg)
    elif optd['cluster_method'] in ['import', 'random', 'skip']:
        pass
    else:
        msg = "Unrecognised cluster_method: {0}".format(optd['cluster_method'])
        exit_util.exit_error(msg)
    if not optd['theseus_exe']:
        optd['theseus_exe'] = 'theseus' + ample_util.EXE_EXT
    try:
        optd['theseus_exe'] = ample_util.find_exe(optd['theseus_exe'])
    except Exception:
        msg = "Cannot find theseus executable: {0}".format(optd['theseus_exe'])
        exit_util.exit_error(msg)
    #
    # SCRWL - we always check for SCRWL as if we are processing QUARK models we want to add sidechains to them
    #
    #if optd['use_scwrl']:
    if not optd['scwrl_exe']:
        optd['scwrl_exe'] = 'Scwrl4' + ample_util.EXE_EXT
    try:
        optd['scwrl_exe'] = ample_util.find_exe(optd['scwrl_exe'])
    except Exception as e:
        LOGGER.info("Cannot find Scwrl executable: {0}".format(optd['scwrl_exe']))
        if optd['use_scwrl']: raise(e)
    #
    # We use shelxe by default so if we can't find it we just warn and set use_shelxe to False
    #
    if optd['use_shelxe']:
        if not optd['shelxe_exe']:
            optd['shelxe_exe'] = 'shelxe' + ample_util.EXE_EXT
        try:
            optd['shelxe_exe'] = ample_util.find_exe(optd['shelxe_exe'])
        except Exception:
            msg = """*** Cannot find shelxe executable in PATH - turning off use of SHELXE. ***
    SHELXE is recommended for the best chance of success. We recommend you install shelxe from:
    http://shelx.uni-ac.gwdg.de/SHELX/
    and install it in your PATH so that AMPLE can use it.
    """
            LOGGER.warn(msg)
            optd['use_shelxe'] = False
    #
    # If shelxe_rebuild is set we need use_shelxe to be set
    #
    if optd['shelxe_rebuild'] and not optd['use_shelxe']:
        msg = 'shelxe_rebuild is set but use_shelxe is False. Please make sure you have shelxe installed.'
        exit_util.exit_error(msg)
    
    if optd['make_frags']:
        if optd['use_homs']:
            LOGGER.info('Making fragments (including homologues)')
        else:
            LOGGER.info('Making fragments EXCLUDING HOMOLOGUES')
    else:
        LOGGER.info('NOT making Fragments')
    
    if optd['make_models']:
        LOGGER.info('\nMaking Rosetta Models')
    else:
        LOGGER.info('NOT making Rosetta Models')
        
        # Print out what is being done
    if optd['use_buccaneer']:
        LOGGER.info('Rebuilding in Bucaneer')
    else:
        LOGGER.info('Not rebuilding in Bucaneer')
    
    if optd['use_arpwarp']:
        LOGGER.info('Rebuilding in ARP/wARP')
    else:
        LOGGER.info('Not rebuilding in ARP/wARP')
    
    # cluster queueing
    if optd['submit_cluster'] and not optd['submit_qtype']:
        msg = 'Must use -submit_qtype argument to specify queueing system (e.g. QSUB, LSF ) if submitting to a cluster.'
        exit_util.exit_error(msg)
    
    if optd['purge']:
        LOGGER.info('*** Purge mode specified - all intermediate files will be deleted ***')
    
    return

def process_restart_options(optd):
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
    if not optd['restart_pkl']: return optd
    if not os.path.isfile(optd['restart_pkl']):
        msg = 'Cannot find restart_pkl file: {0}'.format(optd['restart_pkl'])
        exit_util.exit_error(msg)
    
    LOGGER.info('Restarting from existing pkl file: {0}'.format(optd['restart_pkl']))
    # We use the old dictionary, but udpate it with any new values
    with open(optd['restart_pkl']) as f: optd_old = cPickle.load(f)
    
    # Update key variables that differ with a new run - everything else uses the old values
    optd_old['ample_log'] = optd['ample_log']
    optd_old['run_dir'] = optd['run_dir']
    optd_old['work_dir'] = optd['work_dir']
    optd_old['benchmark_mode'] = optd['benchmark_mode']
    optd_old['benchmark_dir'] = os.path.join(optd['work_dir'], "benchmark")
    optd_old['results_path'] = os.path.join(optd['work_dir'],'resultsd.pkl')
    
    # Now update any variables that were given on the command-line
    for k in optd['cmdline_flags']:
        LOGGER.debug("Restart updating amopt variable: {0} : {1}".format(k, optd[k]))
        optd_old[k] = optd[k]
    
    # We can now replace the old dictionary with this new one
    optd = optd_old
    
    # Go through and see what we need to do
    
    # Reset all variables for doing stuff - otherwise we will always restart from the earliest point
    optd['make_ensembles'] = False
    optd['import_ensembles'] = False # Needs thinking about - have to set so we don't just reimport models/ensembles
    optd['import_models'] = False # Needs thinking about
    optd['make_models'] = False
    optd['make_frags'] = False
    
    # First see if we should benchmark this job. The user may not have supplied a native_pdb with the original
    # job and we only set benchmark mode on seeing the native_pdb
    if optd['native_pdb']:
        if not os.path.isfile(optd['native_pdb']):
            msg = "Cannot find native_pdb: {0}".format(optd['native_pdb'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
        optd['benchmark_mode'] = True
        LOGGER.info('Restart using benchmark mode')
        
    # We always check first to see if there are any mrbump jobs
    optd['mrbump_scripts'] = []
    if 'mrbump_dir' in optd:
        optd['mrbump_scripts'] = mrbump_util.unfinished_scripts(optd)
        if not optd['mrbump_scripts']:
            optd['do_mr'] = False

    if optd['do_mr']:
        if len(optd['mrbump_scripts']):
            LOGGER.info('Restarting from unfinished mrbump scripts: {0}'.format(optd['mrbump_scripts']))
            # Purge unfinished jobs
            for spath in optd['mrbump_scripts']:
                directory, script = os.path.split(spath)
                name, _ = os.path.splitext(script)
                # Hack to delete old job directories
                logfile = os.path.join(directory,name + '.log')
                if os.path.isfile(logfile): os.unlink(logfile)
                jobdir = os.path.join(directory,'search_' + name + '_mrbump')
                if os.path.isdir(jobdir): shutil.rmtree(jobdir)
        elif 'ensembles' in optd and optd['ensembles'] and len(optd['ensembles']):
            # Rerun from ensembles - check for data/ensembles are ok?
            LOGGER.info('Restarting from existing ensembles: {0}'.format(optd['ensembles']))
        elif optd['models_dir'] and optd['models_dir'] and os.path.isdir(optd['models_dir']):
            LOGGER.info('Restarting from existing models: {0}'.format(optd['models_dir']))
            # Check the models
            allsame = False if optd['homologs'] else True 
            if not pdb_edit.check_pdb_directory(optd['models_dir'], sequence=None, single=True, allsame=allsame):
                msg = "Error importing restart models: {0}".format(optd['models_dir'])
                exit_util.exit_error(msg)
            optd['make_ensembles'] = True
        elif optd['frags_3mers'] and optd['frags_9mers']:
            LOGGER.info('Restarting from existing fragments: {0}, {1}'.format(optd['frags_3mers'], optd['frags_9mers']))
            optd['make_models'] = True
    
    return optd

def process_rosetta_options(optd):
    # Create the rosetta modeller - this runs all the checks required
    rosetta_modeller = None
    if optd['make_models'] or optd['make_frags']:  # only need Rosetta if making models
        LOGGER.info('Using ROSETTA so checking options')
        try:
            rosetta_modeller = rosetta_model.RosettaModel(optd=optd)
        except Exception, e:
            msg = "Error setting ROSETTA options: {0}".format(e)
            exit_util.exit_error(msg)
    return rosetta_modeller
