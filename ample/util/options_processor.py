"""Module coordinating the option checking"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Nov 2016"
__version__ = "1.0"

import glob
import logging
import os
import shutil
import sys

from ample.constants import AMPLE_PKL
from ample.ensembler.constants import  SUBCLUSTER_RADIUS_THRESHOLDS, SIDE_CHAIN_TREATMENTS, \
    ALLOWED_SIDE_CHAIN_TREATMENTS, SPICKER_RMSD, SPICKER_TM, POLYALA, RELIABLE, ALLATOM
from ample.modelling import rosetta_model
from ample.util import ample_util
from ample.util import contact_util
from ample.util import exit_util
from ample.util import maxcluster
from ample.util import mrbump_util
from ample.util import mtz_util
from ample.util import pdb_edit
from ample.util import sequence_util

logger = logging.getLogger(__name__)


def check_mandatory_options(optd):
    """Check the mandatory options for correctness

    Description
    -----------
    We check there here rather then with argparse as there doesn't seem
    to be an easy way to get the logic to work of having overlapping
    required and mutually exclusive options

    """

    def _exit(msg, wdir):
        raise RuntimeError(msg)

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

def process_benchmark_options(optd):
    # Benchmark Mode
    if optd['native_pdb'] or optd['benchmark_mode']:
        if optd['native_pdb'] and not os.path.isfile(optd['native_pdb']):
            raise RuntimeError("Cannot find crystal structure PDB: {0}".format(optd['native_pdb']))
        optd['benchmark_mode'] = True
        optd['benchmark_dir'] = os.path.join(optd['work_dir'], "benchmark")
        logger.info("*** AMPLE running in benchmark mode ***")
        # See if we can find TMscore
        if not optd['tmscore_exe']:
            optd['tmscore_exe'] = 'TMscore' + ample_util.EXE_EXT
        try:
            optd['tmscore_exe'] = ample_util.find_exe(optd['tmscore_exe'])
        except ample_util.FileNotFoundError:
            logger.warning("Cannot find TMScore executable: %s", optd['tmscore_exe'])
            optd['maxcluster_exe'] = maxcluster.find_maxcluster(optd)
            optd['have_tmscore'] = False
        else:
            optd['have_tmscore'] = True
            
def process_ensemble_options(optd):
    from ample.ensembler.truncation_util import TRUNCATION_METHODS
    if optd['single_model_mode'] and optd['truncation_scorefile'] and optd['truncation_scorefile_header']:
        #optd['truncation_method'] = "scores"
        optd['truncation_method'] = TRUNCATION_METHODS.SCORES
    elif optd['percent_fixed_intervals']:
        optd['truncation_method'] = TRUNCATION_METHODS.PERCENT_FIXED
    else:
        try:
            optd['truncation_method'] = TRUNCATION_METHODS(optd['truncation_method'])
        except ValueError:
            raise RuntimeError("{} is not a valid truncation method. Use one of: {}".format(optd['truncation_method'],
                                                                                            [e.value for e in TRUNCATION_METHODS]))
        
    # Check we can find all the required programs
    # Maxcluster handled differently as we may need to download the binary
    if optd['subcluster_program'] == 'maxcluster':
        optd['maxcluster_exe'] = maxcluster.find_maxcluster(optd)
    elif optd['subcluster_program'] == 'gesamt':
        if not optd['gesamt_exe']:
            optd['gesamt_exe'] = os.path.join(os.environ['CCP4'], 'bin', 'gesamt' + ample_util.EXE_EXT)
        try:
            optd['gesamt_exe'] = ample_util.find_exe(optd['gesamt_exe'])
        except ample_util.FileNotFoundError as e:
            raise RuntimeError("Cannot find Gesamt executable: {0}".format(optd['gesamt_exe']))
    # Ensemble options
    if optd['cluster_method'] in [SPICKER_RMSD, SPICKER_TM]:
        if not optd['spicker_exe']:
            if optd['cluster_method'] == SPICKER_TM and optd['nproc'] > 1:
                # We need to use the multicore version of SPICKER
                optd['spicker_exe'] = 'spicker_omp' + ample_util.EXE_EXT
            else:
                optd['spicker_exe'] = 'spicker' + ample_util.EXE_EXT
        try:
            optd['spicker_exe'] = ample_util.find_exe(optd['spicker_exe'])
        except ample_util.FileNotFoundError:
            msg = "Cannot find spicker executable: {0}".format(optd['spicker_exe'])
            exit_util.exit_error(msg)
    elif optd['cluster_method'] in ['fast_protein_cluster']:
        if not optd['fast_protein_cluster_exe']:
            optd['fast_protein_cluster_exe'] = 'fast_protein_cluster'
        try:
            optd['fast_protein_cluster_exe'] = ample_util.find_exe(optd['fast_protein_cluster_exe'])
        except ample_util.FileNotFoundError:
            msg = "Cannot find fast_protein_cluster executable: {0}".format(optd['fast_protein_cluster_exe'])
            exit_util.exit_error(msg)
    elif optd['cluster_method'] in ['import', 'random', 'skip']:
        pass
    else:
        raise RuntimeError("Unrecognised cluster_method: {0}".format(optd['cluster_method']))
    if not optd['theseus_exe']:
        optd['theseus_exe'] = os.path.join(os.environ['CCP4'], 'bin', 'theseus' + ample_util.EXE_EXT)
    try:
        optd['theseus_exe'] = ample_util.find_exe(optd['theseus_exe'])
    except ample_util.FileNotFoundError:
        raise RuntimeError("Cannot find theseus executable: {0}".format(optd['theseus_exe']))
    # SCRWL - we always check for SCRWL as if we are processing QUARK models we want to add sidechains to them
    if not optd['scwrl_exe']:
        optd['scwrl_exe'] = os.path.join(os.environ['CCP4'], 'bin', 'Scwrl4' + ample_util.EXE_EXT)
    try:
        optd['scwrl_exe'] = ample_util.find_exe(optd['scwrl_exe'])
    except ample_util.FileNotFoundError as e:
        logger.info("Cannot find Scwrl executable: %s", optd['scwrl_exe'])
        if optd['use_scwrl']:
            raise (e)
    if "subcluster_radius_thresholds" in optd and not optd["subcluster_radius_thresholds"]:
        optd["subcluster_radius_thresholds"] = SUBCLUSTER_RADIUS_THRESHOLDS
    # REM: This should really be disentangled and moved up to definition of all homologs options
    # REM: but could cause confusion with defaults down here.
    if "side_chain_treatments" in optd and not optd["side_chain_treatments"]:
        if optd["homologs"]:
            optd["side_chain_treatments"] = [POLYALA, RELIABLE, ALLATOM]
        else:
            optd["side_chain_treatments"] = SIDE_CHAIN_TREATMENTS
    else:
        optd["side_chain_treatments"] = map(str.lower, optd["side_chain_treatments"])
    unrecognised_sidechains = set(optd["side_chain_treatments"]) - set(ALLOWED_SIDE_CHAIN_TREATMENTS)
    if unrecognised_sidechains:
        raise("Unrecognised side_chain_treatments: {0}".format(unrecognised_sidechains))
    return


def process_options(optd):
    """Process the initial options from the command-line/ample.ini file to set any additional options.

    Description
    -----------
    This is where we take the options determining the type of run we are undertaking and set any additional
    options required based on that runtype. All the major
    """
    # Path for pickling results
    optd['results_path'] = os.path.join(optd['work_dir'], AMPLE_PKL)
    # FASTA processing
    # Check to see if mr_sequence was given and if not mr_sequence defaults to fasta
    if optd['mr_sequence'] != None:
        if not (os.path.exists(str(optd['mr_sequence']))):
            raise RuntimeError('Cannot find mr sequence file: {0}'.format(optd['mr_sequence']))
    else:
        optd['mr_sequence'] = optd['fasta']
    # Process the fasta file and run all the checks on the sequence
    sequence_util.process_fasta(optd, canonicalise=True)
    
    # Not sure if name actually required - see make_fragments.pl
    if optd['name'] and len(optd['name']) != 4:
        raise RuntimeError('-name argument is the wrong length, use 4 chars eg ABCD')
    # Underscore required by rosetta make_fragments.pl
    optd['name'] += '_'
    # MTZ file processing
    try:
        mtz_util.processReflectionFile(optd)
    except Exception, e:
        raise RuntimeError("Error processing reflection file: {0}".format(e))
    # Contact file processing
    if optd['contact_file'] or optd['bbcontacts_file'] or not optd["no_contact_prediction"]:
        contact_util.ContactUtil.check_options(optd)
        optd['use_contacts'] = True
    process_modelling_options(optd)
    # Misc options
    if optd['missing_domain']:
        logger.info('Processing missing domain\n')
        if not os.path.exists(optd['domain_all_chains_pdb']):
            raise RuntimeError('Cannot find file domain_all_chains_pdb: {0}'.format(optd['domain_all_chains_pdb']))
    process_ensemble_options(optd)
    process_mr_options(optd)
    process_benchmark_options(optd)

    logger.info('Running on %d processors', optd['nproc'])
    # cluster queueing
    if optd['submit_qtype']:
        optd['submit_qtype'] = optd['submit_qtype'].upper()
    if optd['submit_cluster'] and not optd['submit_qtype']:
        raise RuntimeError('Must use -submit_qtype argument to specify queueing system (e.g. QSUB, LSF ) if submitting to a cluster.')
    try:
        optd['purge'] = int(optd['purge'])
    except (ValueError, KeyError):
        raise RuntimeError('Purge must be specified as an integer, got: {}'.format(optd['purge']))
    if optd['purge'] > 0:
        logger.info('*** Purge mode level %d specified - intermediate files will be deleted ***', optd['purge'])
    return

def process_modelling_options(optd):
    """ Modelling and ensemble options"""
    # Set default name for modelling directory
    optd['models_dir'] = os.path.join(optd['work_dir'], "models")
    # Check if importing ensembles
    if optd['ensembles']:
        # checks are made in ensembles.import_ensembles
        optd['import_ensembles'] = True
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['cluster_dir']:
        if not os.path.isdir(optd['cluster_dir']):
            raise RuntimeError("Import cluster cannot find directory: {0}".format(optd['cluster_dir']))
        if not glob.glob(os.path.join(optd['cluster_dir'], "*.pdb")):
            raise RuntimeError("Import cluster cannot find pdbs in directory: {0}".format(optd['cluster_dir']))
        logger.info("Importing pre-clustered models from directory: %d\n", optd['cluster_dir'])
        optd['cluster_method'] = 'import'
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['ideal_helices']:
        optd['make_frags'] = False
        optd['make_models'] = False
    elif optd['homologs']:
        optd['make_frags'] = False
        optd['make_models'] = False
        if not os.path.isfile(str(optd['alignment_file'])):
            # We need to use gesamt or mustang to do the alignment
            if optd['homolog_aligner'] == 'gesamt':
                if not ample_util.is_exe(str(optd['gesamt_exe'])):
                    optd['gesamt_exe'] = os.path.join(os.environ['CCP4'], 'bin', 'gesamt' + ample_util.EXE_EXT)
                if not ample_util.is_exe(str(optd['gesamt_exe'])):
                    raise RuntimeError('Using homologs without an alignment file and cannot find gesamt_exe: {0}'.format(
                        optd['gesamt_exe']))
            elif optd['homolog_aligner'] == 'mustang':
                if not ample_util.is_exe(str(optd['mustang_exe'])):
                    raise RuntimeError('Using homologs without an alignment file and cannot find mustang_exe: {0}'.format(
                        optd['mustang_exe']))
            else:
                raise RuntimeError('Unknown homolog_aligner: {0}'.format(optd['homolog_aligner']))
        if not os.path.isdir(str(optd['models'])):
            raise RuntimeError("Homologs option requires a directory of pdb models to be supplied\n" + \
                              "Please supply the models with the -models flag")
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
    # Check import flags
    if optd['import_ensembles'] and (optd['import_models']):
        raise RuntimeError("Cannot import both models and ensembles/clusters!")
    # NMR Checks
    if optd['nmr_model_in']:
        logger.info("Using nmr_model_in file: %s", optd['nmr_model_in'])
        if not os.path.isfile(optd['nmr_model_in']):
            msg = "nmr_model_in flag given, but cannot find file: {0}".format(optd['nmr_model_in'])
            exit_util.exit_error(msg)
        if optd['nmr_remodel']:
            optd['make_models'] = True
            if optd['nmr_remodel_fasta']:
                if not os.path.isfile(optd['nmr_remodel_fasta']):
                    raise RuntimeError("Cannot find nmr_remodel_fasta file: {0}".format(optd['nmr_remodel_fasta']))
            else:
                optd['nmr_remodel_fasta'] = optd['fasta']
            msg = "NMR model will be remodelled with ROSETTA using the sequence from: {0}".format(
                optd['nmr_remodel_fasta'])
            logger.info(msg)
            if not (optd['frags_3mers'] and optd['frags_9mers']):
                optd['make_frags'] = True
                msg = "nmr_remodel - will be making our own fragment files"
                logger.info(msg)
            else:
                if not (os.path.isfile(optd['frags_3mers']) and os.path.isfile(optd['frags_9mers'])):
                    raise RuntimeError("frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(
                        optd['frags_3mers'], optd['frags_9mers']))
                optd['make_frags'] = False
        else:
            optd['make_frags'] = False
            optd['make_models'] = False
            msg = "Running in NMR truncate only mode"
            logger.info(msg)
    elif optd['make_models']:
        if not os.path.isdir(optd['models_dir']):
            os.mkdir(optd['models_dir'])
        # If the user has given both fragment files we check they are ok and unset make_frags
        if optd['frags_3mers'] and optd['frags_9mers']:
            if not os.path.isfile(optd['frags_3mers']) or not os.path.isfile(optd['frags_9mers']):
                raise RuntimeError("frags_3mers and frag_9mers files given, but cannot locate them:\n{0}\n{1}\n".format(
                    optd['frags_3mers'], optd['frags_9mers']))
            optd['make_frags'] = False
        if optd['make_frags'] and (optd['frags_3mers'] or optd['frags_9mers']):
            raise RuntimeError("make_frags set to true, but you have given the path to the frags_3mers or frags_9mers")
        if not optd['make_frags'] and not (optd['frags_3mers'] and optd['frags_9mers']):
            msg = """*** Missing fragment files! ***
    Please supply the paths to the fragment files using the -frags_3mers and -frags_9mers flags.
    These can be generated using the Robetta server: http://robetta.bakerlab.org
    Please see the AMPLE documentation for further information."""
            raise RuntimeError(msg)
    if optd['make_frags']:
        if optd['use_homs']:
            logger.info('Making fragments (including homologues)')
        else:
            logger.info('Making fragments EXCLUDING HOMOLOGUES')
    else:
        logger.info('NOT making Fragments')
    if optd['make_models']:
        logger.info('\nMaking Rosetta Models')
    else:
        logger.info('NOT making Rosetta Models')


def process_mr_options(optd):
    # Molecular Replacement Options
    if optd['molrep_only']:
        optd['phaser_only'] = False
    if optd['molrep_only']:
        optd['mrbump_programs'] = ['molrep']
    elif optd['phaser_only']:
        optd['mrbump_programs'] = ['phaser']
    else:
        optd['mrbump_programs'] = ['molrep', 'phaser']
    if optd['phaser_rms'] != 'auto':
        try:
            phaser_rms = float(optd['phaser_rms'])
        except ValueError as e:
            msg = "Error converting phaser_rms '{0}' to floating point: {1}".format(optd['phaser_rms'], e)
            exit_util.exit_error(msg)
        else:
            optd['phaser_rms'] = phaser_rms
    # We use shelxe by default so if we can't find it we just warn and set use_shelxe to False
    if optd['use_shelxe']:
        if optd['mtz_min_resolution'] > mrbump_util.SHELXE_MAX_PERMITTED_RESOLUTION:
            logger.warn("Disabling use of SHELXE as min resolution of %f is < accepted limit of %f",
                            optd['mtz_min_resolution'],
                            mrbump_util.SHELXE_MAX_PERMITTED_RESOLUTION)
            optd['use_shelxe'] = False
            optd['shelxe_rebuild'] = False
    if optd['use_shelxe']:
        if not optd['shelxe_exe']:
            optd['shelxe_exe'] = os.path.join(os.environ['CCP4'], 'bin', 'shelxe' + ample_util.EXE_EXT)
        try:
            optd['shelxe_exe'] = ample_util.find_exe(optd['shelxe_exe'])
        except ample_util.FileNotFoundError:
            msg = """*** Cannot find shelxe executable in PATH - turning off use of SHELXE. ***
    SHELXE is recommended for the best chance of success. We recommend you install shelxe from:
    http://shelx.uni-ac.gwdg.de/SHELX/
    and install it in your PATH so that AMPLE can use it.
    """
            logger.warn(msg)
            optd['use_shelxe'] = False
    if optd['shelxe_rebuild']:
        optd['shelxe_rebuild_arpwarp'] = True
        optd['shelxe_rebuild_buccaneer'] = True

    if optd['refine_rebuild_arpwarp'] or optd['shelxe_rebuild_arpwarp']:
        auto_tracing_sh = None
        if 'warpbin' in os.environ:
            _path = os.path.join(os.environ['warpbin'], "auto_tracing.sh")
            if os.path.isfile(_path):
                auto_tracing_sh = _path
        if auto_tracing_sh:
            logger.info('Using arpwarp script: %s', auto_tracing_sh)
        else:
            logger.warn('Cannot find arpwarp script! Disabling use of arpwarp.')
            optd['refine_rebuild_arpwarp'] = False
            optd['shelxe_rebuild_arpwarp'] = False

    if optd['refine_rebuild_arpwarp'] or optd['shelxe_rebuild_arpwarp']:
        logger.info('Rebuilding in ARP/wARP')
    else:
        logger.info('Not rebuilding in ARP/wARP')

    if optd['refine_rebuild_buccaneer'] or optd['shelxe_rebuild_buccaneer']:
        logger.info('Rebuilding in Bucaneer')
    else:
        logger.info('Not rebuilding in Bucaneer')

    # If shelxe_rebuild is set we need use_shelxe to be set
    if optd['shelxe_rebuild'] and not optd['use_shelxe']:
        raise RuntimeError('shelxe_rebuild is set but use_shelxe is False. Please make sure you have shelxe installed.')


def process_restart_options(optd):
    """Process the restart options

    Description
    -----------
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

    Notes
    -----
    We return the dictionary as we may need to change it and it seems we can't change the external
    reference in this scope. I think?...

    """
    if not optd['restart_pkl']:
        return optd
    logger.info('Restarting from existing pkl file: %s', optd['restart_pkl'])

    # Go through and see what we need to do
    # Reset all variables for doing stuff - otherwise we will always restart from the earliest point
    optd['make_ensembles'] = False
    #optd['import_ensembles'] = False # Needs thinking about - have to set so we don't just reimport models/ensembles
    optd['import_models'] = False  # Needs thinking about
    optd['make_models'] = False
    optd['make_frags'] = False

    # First see if we should benchmark this job. The user may not have supplied a native_pdb with the original
    # job and we only set benchmark mode on seeing the native_pdb
    if optd['native_pdb']:
        if not os.path.isfile(optd['native_pdb']):
            raise RuntimeError("Cannot find native_pdb: {0}".format(optd['native_pdb']))
        optd['benchmark_mode'] = True
        logger.info('Restart using benchmark mode')

    # We always check first to see if there are any mrbump jobs
    optd['mrbump_scripts'] = []
    if 'mrbump_dir' in optd:
        optd['mrbump_scripts'] = mrbump_util.unfinished_scripts(optd)
        if not optd['mrbump_scripts']:
            optd['do_mr'] = False

    if optd['do_mr']:
        if len(optd['mrbump_scripts']):
            logger.info('Restarting from unfinished mrbump scripts: %s', optd['mrbump_scripts'])
            # Purge unfinished jobs
            for spath in optd['mrbump_scripts']:
                directory, script = os.path.split(spath)
                name, _ = os.path.splitext(script)
                # Hack to delete old job directories
                logfile = os.path.join(directory, name + '.log')
                if os.path.isfile(logfile):
                    os.unlink(logfile)
                jobdir = os.path.join(directory, 'search_' + name + '_mrbump')
                if os.path.isdir(jobdir):
                    shutil.rmtree(jobdir)
        elif 'ensembles' in optd and optd['ensembles'] and len(optd['ensembles']):
            # Rerun from ensembles - check for data/ensembles are ok?
            logger.info('Restarting from existing ensembles: %s', optd['ensembles'])
        elif 'models_dir' in optd and optd['models_dir'] and os.path.isdir(optd['models_dir']):
            logger.info('Restarting from existing models: %s', optd['models_dir'])
            allsame = False if optd['homologs'] else True
            if not pdb_edit.check_pdb_directory(optd['models_dir'], sequence=None, single=True, allsame=allsame):
                raise RuntimeError("Error importing restart models: {0}".format(optd['models_dir']))
            optd['make_ensembles'] = True
        elif optd['frags_3mers'] and optd['frags_9mers']:
            logger.info('Restarting from existing fragments: %s, %s', optd['frags_3mers'], optd['frags_9mers'])
            optd['make_models'] = True
    return optd


def process_rosetta_options(optd):
    # Create the rosetta modeller - this runs all the checks required
    rosetta_modeller = None
    if optd['make_models'] or optd['make_frags']:  # only need Rosetta if making models
        logger.info('Using ROSETTA so checking options')
        try:
            rosetta_modeller = rosetta_model.RosettaModel(optd=optd)
        except Exception as e:
            msg = "Error setting ROSETTA options: {0}".format(e)
            exit_util.exit_error(msg)
        optd['modelling_workdir'] = rosetta_modeller.work_dir
    return rosetta_modeller


def restart_amoptd(optd):
    """Create an ample dictionary from a restart pkl file

    Description
    -----------
    For any new command-line options, we update the old dictionary with the new values
    We then go through the new dictionary and set ant of the flags corresponding to the data we find:

    Notes
    -----
    We return the dictionary as we may need to change it and it seems we can't change the external
    reference in this scope. I think?...

    """
    if not optd['restart_pkl']:
        return optd
    logger.info('Restarting from existing pkl file: %s', optd['restart_pkl'])
    optd_old = ample_util.read_amoptd(optd['restart_pkl'])
    for k in optd['cmdline_flags']:
        logger.debug("Restart updating amopt variable: %s : %s", k, str(optd[k]))
        optd_old[k] = optd[k]
    optd = optd_old
    return optd
