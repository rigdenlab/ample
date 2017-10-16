#!/usr/bin/env ccp4-python
"""
This is AMPLE
"""
import argparse
import logging
import os
import platform
import shutil
import sys
import time

from ample import ensembler
from ample.ensembler.constants import UNMODIFIED
from ample.util import ample_util
from ample.util import argparse_util
from ample.util import benchmark_util
from ample.util import config_util
from ample.util import contact_util
from ample.util import exit_util
from ample.util import logging_util
from ample.util import mrbump_util
from ample.util import options_processor
from ample.util import pdb_edit
from ample.util import pyrvapi_results
from ample.util import workers_util
from ample.util import version

__author__ = "Jens Thomas, Felix Simkovic, Adam Simpkin, Ronan Keegan, and Jaclyn Bibby"
__credits__ = "Daniel Rigden, Martyn Winn, and Olga Mayans"
__email__ = "drigden@liverpool.ac.uk"
__version__ = version.__version__

logger = logging_util.setup_console_logging()
monitor = None

class Ample(object):
    """Class to generate ensembles from ab inito models (all models must have same sequence)
    
    """

    def __init__(self):
        self.amopt = None
        self.ample_output = None
        return

    def main(self, args=None):
        """Main AMPLE routine.
        
        We require this as the multiprocessing module (only on **!!*%$$!! Windoze) 
        requires that the main module can be imported. We there need ample to be 
        a python script that can be imported, hence the main routine with its 
        calling protected by the if __name__=="__main__":...
        
        args is an option argument that can contain the command-line arguments 
        for the program - required for testing.
        """
        argso = self.process_command_line(args=args)

        self.amopt = amopt = config_util.AMPLEConfigOptions()
        amopt.populate(argso)

        # Setup things like logging, file structure, etc...
        amopt.d = self.setup(amopt.d)
        rosetta_modeller = options_processor.process_rosetta_options(amopt.d)

        # Display the parameters used
        logger.debug(amopt.prettify_parameters())

        amopt.write_config_file()
        #######################################################
        # SCRIPT PROPER STARTS HERE
        time_start = time.time()

        # Create function for monitoring jobs - static function decorator?
        if self.ample_output:
            def monitor():
                return self.ample_output.display_results(amopt.d)
        else:
            monitor = None

        if amopt.d['benchmark_mode'] and amopt.d['native_pdb']:
            # Process the native before we do anything else
            benchmark_util.analysePdb(amopt.d)
            
        # Create constituent models from an NMR ensemble
        if amopt.d['nmr_model_in']:
            logger.info('Splitting NMR ensemble into constituent models')
            amopt.d['models'] = pdb_edit.split_pdb(amopt.d['nmr_model_in'],
                                                   directory=amopt.d['models_dir'],
                                                   strip_hetatm=True,
                                                   same_size=True)
            logger.info('NMR ensemble contained {0} models'.format(len(amopt.d['models'])))

        # Modelling business happens here
        if amopt.d['make_frags'] or amopt.d['make_models']:
            self.modelling(amopt.d, rosetta_modeller)
            amopt.write_config_file()

        # Ensembling business next
        if amopt.d['make_ensembles']:
            self.ensembling(amopt.d)
            amopt.write_config_file()

        # Some MR here
        if amopt.d['do_mr']:
            self.molecular_replacement(amopt.d)
            amopt.write_config_file()

        # Timing data
        time_stop = time.time()
        elapsed_time = time_stop - time_start
        run_in_min = elapsed_time / 60
        run_in_hours = run_in_min / 60
        msg = os.linesep + \
            'All processing completed  (in {0:6.2F} hours)'.format(
                run_in_hours) + os.linesep
        msg += '----------------------------------------' + os.linesep
        logging.info(msg)

        # Benchmark mode
        if amopt.d['benchmark_mode']:
            self.benchmarking(amopt.d)
            amopt.write_config_file()

        amopt.write_config_file()
        # Flag to show that we reached the end without error - useful for integration testing
        amopt.d['AMPLE_finished'] = True
        ample_util.save_amoptd(amopt.d)

        logger.info("AMPLE finished at: {0}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        logger.info(ample_util.reference.format(
            refs=ample_util.construct_references(amopt.d)))
        logger.info(ample_util.footer)

        # Finally update pyrvapi results
        if self.ample_output:
            self.ample_output.display_results(amopt.d)
        if self.ample_output:
            self.ample_output.rvapi_shutdown(amopt.d)
        return

    def benchmarking(self, optd):
        if optd['submit_cluster']:
            # Pickle dictionary so it can be opened by the job to get the parameters
            ample_util.save_amoptd(optd)
            script = benchmark_util.cluster_script(optd)
            workers_util.run_scripts(job_scripts=[script],
                                     monitor=monitor,
                                     nproc=optd['nproc'],
                                     job_time=43200,
                                     job_name='benchmark',
                                     submit_cluster=optd['submit_cluster'],
                                     submit_qtype=optd['submit_qtype'],
                                     submit_queue=optd['submit_queue'],
                                     submit_array=optd['submit_array'],
                                     submit_max_array=optd['submit_max_array'])
            # queue finished so unpickle results
            optd.update(ample_util.read_amoptd(optd['results_path']))
        else:
            benchmark_util.analyse(optd)
            ample_util.save_amoptd(optd)

        return

    def ensembling(self, optd):

        if optd['import_ensembles']:
            ensembler.import_ensembles(optd)
        elif optd['ideal_helices']:
            ample_util.ideal_helices(optd)
            logger.info("*** Using ideal helices to solve structure ***")
        else:
            # Import the models here instead of cluster_util.
            if optd['cluster_method'] is 'import':
                # HACK - this is certainly not how we want to do it. One flag for all (-models) in future
                optd['models'] = optd['cluster_dir']
                optd['models'] = ample_util.extract_models(optd)

            # Check we have some models to work with
            if not (optd['single_model_mode'] or optd['models']):
                ample_util.save_amoptd(optd)
                msg = "ERROR! Cannot find any pdb files in: {0}".format(
                    optd['models_dir'])
                exit_util.exit_error(msg)
            optd['ensemble_ok'] = os.path.join(optd['work_dir'], 'ensemble.ok')
            if optd['submit_cluster']:
                # Pickle dictionary so it can be opened by the job to get the parameters
                ample_util.save_amoptd(optd)
                script = ensembler.cluster_script(optd)
                workers_util.run_scripts(job_scripts=[script],
                                         monitor=monitor,
                                         nproc=optd['nproc'],
                                         job_time=optd['ensembler_timeout'],
                                         job_name='ensemble',
                                         submit_cluster=optd['submit_cluster'],
                                         submit_qtype=optd['submit_qtype'],
                                         submit_queue=optd['submit_queue'],
                                         submit_pe_lsf=optd['submit_pe_lsf'],
                                         submit_pe_sge=optd['submit_pe_sge'],
                                         submit_array=optd['submit_array'],
                                         submit_max_array=optd['submit_max_array'])
                # queue finished so unpickle results
                optd.update(ample_util.read_amoptd(optd['results_path']))
            else:
                try:
                    ensembler.create_ensembles(optd)
                except Exception as e:
                    msg = "Error creating ensembles: {0}".format(e)
                    exit_util.exit_error(msg, sys.exc_info()[2])

            # Check we have something to work with
            if not os.path.isfile(optd['ensemble_ok']) or 'ensembles' not in optd.keys() or not len(optd['ensembles']):
                msg = "Problem generating ensembles!"
                exit_util.exit_error(msg)

            if not (optd['homologs'] or optd['single_model_mode']):
                ensemble_summary = ensembler.ensemble_summary(
                    optd['ensembles_data'])
                logger.info(ensemble_summary)

        # Save the results
        ample_util.save_amoptd(optd)

        # Bail here if we didn't create anything
        if not len(optd['ensembles']):
            msg = "### AMPLE FAILED TO GENERATE ANY ENSEMBLES! ###\nExiting..."
            exit_util.exit_error(msg)

        # Update results view
        if self.ample_output:
            self.ample_output.display_results(optd)
        return

    def modelling(self, optd, rosetta_modeller=None):
        if not rosetta_modeller:
            rosetta_modeller = options_processor.process_rosetta_options(optd)
        # Make Rosetta fragments
        if optd['make_frags']:
            rosetta_modeller.generate_fragments(optd)
            optd['frags_3mers'] = rosetta_modeller.frags_3mers
            optd['frags_9mers'] = rosetta_modeller.frags_9mers
            optd['psipred_ss2'] = rosetta_modeller.psipred_ss2

        if optd["use_contacts"] and not optd['restraints_file']:
            con_util = contact_util.ContactUtil(
                optd['fasta'], 'fasta',
                contact_file=optd['contact_file'],
                contact_format=optd['contact_format'],
                bbcontacts_file=optd['bbcontacts_file'],
                bbcontacts_format=optd["bbcontacts_format"],
                cutoff_factor=optd['restraints_factor'],
                distance_to_neighbor=optd['distance_to_neighbour']
            )

            optd["contacts_dir"] = os.path.join(optd["work_dir"], "contacts")
            if not os.path.isdir(optd["contacts_dir"]):
                os.mkdir(optd["contacts_dir"])
            if con_util.require_contact_prediction:
                if con_util.found_ccmpred_contact_prediction_deps:
                    con_util.predict_contacts_from_sequence(
                        wdir=optd["contacts_dir"])
                    optd["contact_file"] = con_util.contact_file
                    optd["contact_format"] = con_util.contact_format

            if con_util.do_contact_analysis:
                plot_file = os.path.join(
                    optd['contacts_dir'], optd['name'] + ".cm.png")
                if optd['native_pdb'] and optd['native_pdb_std']:
                    structure_file = optd['native_pdb_std']
                elif optd["native_pdb"]:
                    structure_file = optd['native_std']
                else:
                    structure_file = None
                optd['contact_map'], optd['contact_ppv'] = con_util.summarize(
                    plot_file, structure_file, 'pdb', optd['native_cutoff'])

                restraints_file = os.path.join(
                    optd['contacts_dir'], optd['name'] + ".cst")
                optd['restraints_file'] = con_util.write_restraints(
                    restraints_file, optd['restraints_format'], optd['energy_function']
                )
            else:
                con_util = None
        else:
            con_util = None

        if optd['make_models'] and optd['restraints_file']:
            rosetta_modeller.restraints_file = optd['restraints_file']

        if optd['make_models']:
            # Make the models
            logger.info('----- making Rosetta models--------')
            if optd['nmr_remodel']:
                try:
                    optd['models'] = rosetta_modeller.nmr_remodel(models=optd['models'],
                                                                  ntimes=optd['nmr_process'],
                                                                  alignment_file=optd['alignment_file'],
                                                                  remodel_fasta=optd['nmr_remodel_fasta'],
                                                                  monitor=monitor)
                except Exception as e:
                    msg = "Error remodelling NMR ensemble: {0}".format(e)
                    exit_util.exit_error(msg, sys.exc_info()[2])
            else:
                logger.info('making {0} models...'.format(optd['nmodels']))
                try:
                    optd['models'] = rosetta_modeller.ab_initio_model(
                        monitor=monitor)
                except Exception as e:
                    msg = "Error running ROSETTA to create models: {0}".format(
                        e)
                    exit_util.exit_error(msg, sys.exc_info()[2])
                if not pdb_edit.check_pdb_directory(optd['models_dir'], sequence=optd['sequence']):
                    msg = "Problem with rosetta pdb files - please check the log for more information"
                    exit_util.exit_error(msg)
                msg = 'Modelling complete - models stored in: {0}\n'.format(
                    optd['models_dir'])
                logger.info(msg)

        elif optd['import_models']:
            logger.info('Importing models from directory: {0}\n'.format(
                optd['models_dir']))
            if optd['homologs']:
                optd['models'] = ample_util.extract_models(
                    optd, sequence=None, single=True, allsame=False)
            else:
                optd['models'] = ample_util.extract_models(optd)
                # Need to check if Quark and handle things accordingly
                if optd['quark_models']:
                    # We always add sidechains to QUARK models if SCWRL is installed
                    if ample_util.is_exe(optd['scwrl_exe']):
                        optd['use_scwrl'] = True
                    else:
                        # No SCWRL so don't do owt with the side chains
                        logger.info('Using QUARK models but SCWRL is not installed '
                                    'so only using {0} sidechains'.format(UNMODIFIED))
                        optd['side_chain_treatments'] = [UNMODIFIED]

        # Sub-select the decoys using contact information
        if con_util and optd['subselect_mode'] and not (optd['nmr_model_in'] or optd['nmr_remodel']):
            logger.info(
                'Subselecting models from directory using provided contact information')
            optd['models'] = con_util.subselect_decoys(
                optd['models'], 'pdb', mode=optd['subselect_mode'], **optd)

        # Save the results
        ample_util.save_amoptd(optd)

        return

    def molecular_replacement(self, optd):

        if not optd['mrbump_scripts']:
            # MRBUMP analysis of the ensembles
            logger.info('----- Running MRBUMP on ensembles--------\n\n')
            if len(optd['ensembles']) < 1:
                msg = "ERROR! Cannot run MRBUMP as there are no ensembles!"
                exit_util.exit_error(msg)

            if optd['mrbump_dir'] is None:
                bump_dir = os.path.join(optd['work_dir'], 'MRBUMP')
                optd['mrbump_dir'] = bump_dir
            else:
                bump_dir = optd['mrbump_dir']
            if not os.path.exists(bump_dir):
                os.mkdir(bump_dir)

            optd['mrbump_results'] = []
            logger.info(
                "Running MRBUMP jobs in directory: {0}".format(bump_dir))

            # Set an ensemble-specific phaser_rms if required
            if optd['phaser_rms'] == 'auto':
                ensembler.set_phaser_rms_from_subcluster_score(optd)

            # Sort the ensembles in a favourable way
            logger.info("Sorting ensembles")
            sort_keys = ['cluster_num', 'truncation_level',
                         'subcluster_radius_threshold', 'side_chain_treatment']
            ensemble_pdbs_sorted = ensembler.sort_ensembles(optd['ensembles'], optd['ensembles_data'],
                                                            keys=sort_keys, prioritise=True)

            # Create job scripts
            logger.info("Generating MRBUMP runscripts")
            optd['mrbump_scripts'] = mrbump_util.write_mrbump_files(ensemble_pdbs_sorted,
                                                                    optd,
                                                                    job_time=mrbump_util.MRBUMP_RUNTIME,
                                                                    ensemble_options=optd['ensemble_options'],
                                                                    directory=bump_dir)

        # Create function for monitoring jobs - static function decorator?
        if self.ample_output:
            def monitor():
                r = mrbump_util.ResultsSummary()
                r.extractResults(optd['mrbump_dir'], purge=optd['purge'])
                optd['mrbump_results'] = r.results
                return self.ample_output.display_results(optd)
        else:
            monitor = None

        # Save results here so that we have the list of scripts and mrbump directory set
        ample_util.save_amoptd(optd)

        # Change to mrbump directory before running
        os.chdir(optd['mrbump_dir'])
        ok = workers_util.run_scripts(job_scripts=optd['mrbump_scripts'],
                                      monitor=monitor,
                                      check_success=mrbump_util.checkSuccess,
                                      early_terminate=optd['early_terminate'],
                                      nproc=optd['nproc'],
                                      job_time=mrbump_util.MRBUMP_RUNTIME,
                                      job_name='mrbump',
                                      submit_cluster=optd['submit_cluster'],
                                      submit_qtype=optd['submit_qtype'],
                                      submit_queue=optd['submit_queue'],
                                      submit_array=optd['submit_array'],
                                      submit_max_array=optd['submit_max_array'])

        if not ok:
            msg = "Error running MRBUMP on the ensembles!\nCheck logs in directory: {0}".format(
                optd['mrbump_dir'])
            exit_util.exit_error(msg)

        # Collect the MRBUMP results
        results_summary = mrbump_util.ResultsSummary()
        optd['mrbump_results'] = results_summary.extractResults(
            optd['mrbump_dir'], purge=optd['purge'])
        optd['success'] = results_summary.success

        ample_util.save_amoptd(optd)

        # Now print out the final summary
        summary = mrbump_util.finalSummary(optd)
        logger.info(summary)

        return

    def process_command_line(self, args=None, contacts=True, modelling=True, mol_rep=True):
        """Process the command-line.
        :args: optional argument that can hold the command-line arguments if we 
        have been called from within python for testing
        """
        parser = argparse.ArgumentParser(description="AMPLE: Ab initio Modelling of Proteins for moLEcular replacement",
                                         prefix_chars="-")
        argparse_util.add_general_options(parser)
        argparse_util.add_cluster_submit_options(parser)
        ensembler.add_argparse_options(parser)

        if contacts:
            argparse_util.add_contact_options(parser)
        if mol_rep:
            argparse_util.add_mr_options(parser)
        if modelling:
            argparse_util.add_rosetta_options(parser)

        return parser.parse_args(args)

    def setup(self, optd):
        """We take and return an ample dictionary as an argument. This is required because options_processor.process_restart_options
        Changes what optd points at, and this means that if we just use the reference, we end up pointing at the old, obselete dictionary"""

        # Update the ample dictionary in case we are restarting
        optd = options_processor.restart_amoptd(optd)

        # Make a work directory - this way all output goes into this directory
        if optd['work_dir'] and not optd['restart_pkl']:
            logger.info(
                'Making a named work directory: {0}'.format(optd['work_dir']))
            try:
                os.mkdir(optd['work_dir'])
            except Exception as e:
                msg = "Cannot create work_dir {0}: {1}".format(
                    optd['work_dir'], e)
                exit_util.exit_error(msg, sys.exc_info()[2])

        if not optd['work_dir']:
            if not os.path.exists(optd['run_dir']):
                msg = 'Cannot find run directory: {0}'.format(optd['run_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])

            if bool(optd['rvapi_document']):
                # With JSCOFE we run in the run directory
                optd['work_dir'] = optd['run_dir']
            else:
                logger.info('Making a run directory: checking for previous runs...')
                optd['work_dir'] = ample_util.make_workdir(optd['run_dir'],
                                                           ccp4i2=bool(optd['ccp4i2_xml']))
        # Go to the work directory
        os.chdir(optd['work_dir'])

        # Set up logging
        ample_log = os.path.join(optd['work_dir'], 'AMPLE.log')
        debug_log = os.path.join(optd['work_dir'], 'debug.log')
        optd['ample_log'] = ample_log

        # Set up ample output file and debug log file.
        logging_util.setup_file_logging(ample_log, level=logging.INFO)
        logging_util.setup_file_logging(debug_log, level=logging.DEBUG)

        # Make sure the CCP4 environment is set up properly
        ccp4_home = self.setup_ccp4(optd)
        ccp4_version = ".".join([str(x) for x in optd['ccp4_version']])

        # Print out Version and invocation
        logger.info(ample_util.header)
        logger.info("AMPLE version: {0}".format(version.__version__))
        logger.info("Running with CCP4 version: {0} from directory: {1}".format(
            ccp4_version, ccp4_home))
        logger.info("Running on host: {0}".format(platform.node()))
        logger.info("Running on platform: {0}".format(platform.platform()))
        logger.info("Job started at: {0}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        logger.info(
            "Invoked with command-line:\n{0}\n".format(" ".join(sys.argv)))
        logger.info("Running in directory: {0}\n".format(optd['work_dir']))

        # Display pyrvapi results
        if pyrvapi_results.pyrvapi:
            self.ample_output = pyrvapi_results.AmpleOutput(optd)
            self.ample_output.display_results(optd)

        # Check mandatory/exclusive options
        options_processor.check_mandatory_options(optd)

        # Check if we are restarting from an existing pkl file - we don't process the options from this
        # run if so
        optd = options_processor.process_restart_options(optd)
        if not optd['restart_pkl']:
            # Only process the remaining options if we aren't in restart mode
            options_processor.process_options(optd)

        # Bail and clean up if we were only checking the options
        if optd['dry_run']:
            logger.info('Dry run finished checking options - cleaning up...')
            os.chdir(optd['run_dir'])
            shutil.rmtree(optd['work_dir'])
            sys.exit(0)

        logger.info('All needed programs are found, continuing...')

        return optd

    def setup_ccp4(self, amoptd):
        """Check CCP4 is available and return the top CCP4 directory"""
        # Make sure CCP4 is around
        if not "CCP4" in os.environ:
            msg = "Cannot find CCP4 installation - please make sure CCP4 is installed and the setup scripts have been run!"
            exit_util.exit_error(msg)

        if not "CCP4_SCR" in os.environ:
            msg = "$CCP4_SCR environement variable not set - please make sure CCP4 is installed and the setup scripts have been run!"
            exit_util.exit_error(msg)

        if not os.path.isdir(os.environ['CCP4_SCR']):
            msg = "*** WARNING ***\n"
            msg += "Cannot find the $CCP4_SCR directory: {0}\n".format(
                os.environ['CCP4_SCR'])
            msg += "The directory will be created, but it should have already been created by the CCP4 startup scripts\n"
            msg += "Please make sure CCP4 is installed and the setup scripts have been run."
            logger.critical(msg)
            os.mkdir(os.environ['CCP4_SCR'])
            #exit_util.exit_error(msg)

        # Record the CCP4 version we're running with  - also required in pyrvapi_results
        amoptd['ccp4_version'] = ample_util.ccp4_version()

        return os.environ['CCP4']


if __name__ == "__main__":
    try:
        Ample().main()
    except Exception as e:
        msg = "Error running main AMPLE program: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
