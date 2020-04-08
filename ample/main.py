#!/usr/bin/env ccp4-python
"""
This is AMPLE
"""
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
from ample.util import process_models
from ample.util import pyrvapi_results
from ample.util import reference_manager
from ample.util import version

from pyjob.factory import TaskFactory

__author__ = "Jens Thomas, Felix Simkovic, Adam Simpkin, Ronan Keegan, and Jaclyn Bibby"
__credits__ = "Daniel Rigden, Martyn Winn, and Olga Mayans"
__email__ = "drigden@liverpool.ac.uk"
__version__ = version.__version__

logger = None
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
        argso = argparse_util.process_command_line(args=args)
        # Work directory and loggers need to be setup before we do anything else
        self.setup_workdir(argso)
        global logger
        logger = logging_util.setup_logging(argso)

        # Logging and work directories in place so can start work
        self.amopt = amopt = config_util.AMPLEConfigOptions()
        amopt.populate(argso)
        amopt.d = self.setup(amopt.d)
        rosetta_modeller = options_processor.process_rosetta_options(amopt.d)
        logger.debug(amopt.prettify_parameters())  # Display the parameters used
        amopt.write_config_file()
        time_start = time.time()
        if self.ample_output:

            def monitor():
                return self.ample_output.display_results(amopt.d)

        else:
            monitor = None

        # Process any files we may have been given
        model_results = process_models.extract_and_validate_models(amopt.d)
        if model_results:
            process_models.handle_model_import(amopt.d, model_results)
        if amopt.d['benchmark_mode'] and amopt.d['native_pdb']:
            # Process the native before we do anything else
            benchmark_util.analysePdb(amopt.d)

        # Create constituent models from an NMR ensemble
        if amopt.d['nmr_model_in']:
            nmr_mdir = os.path.join(amopt.d['work_dir'], 'nmr_models')
            amopt.d['modelling_workdir'] = nmr_mdir
            logger.info('Splitting NMR ensemble into constituent models in directory: {0}'.format(nmr_mdir))
            amopt.d['processed_models'] = pdb_edit.split_pdb(
                amopt.d['nmr_model_in'], directory=nmr_mdir, strip_hetatm=True, same_size=True
            )
            logger.info('NMR ensemble contained {0} models'.format(len(amopt.d['processed_models'])))

        # Modelling business happens here
        if self.modelling_required(amopt.d):
            self.modelling(amopt.d, rosetta_modeller)
            ample_util.save_amoptd(amopt.d)
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
        msg = os.linesep + 'All processing completed  (in {0:6.2F} hours)'.format(run_in_hours) + os.linesep
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

        logger.info("AMPLE finished at: %s", time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
        ref_mgr = reference_manager.ReferenceManager(amopt.d)
        ref_mgr.save_citations_to_file(amopt.d)
        logger.info(ref_mgr.citations_as_text)
        logger.info(reference_manager.footer)

        # Finally update pyrvapi results
        if self.ample_output:
            self.ample_output.display_results(amopt.d)
            self.ample_output.rvapi_shutdown(amopt.d)

        self.cleanup(amopt.d)
        return

    def benchmarking(self, optd):
        if optd['submit_cluster']:
            # Pickle dictionary so it can be opened by the job to get the parameters
            ample_util.save_amoptd(optd)
            script = benchmark_util.cluster_script(optd)
            with TaskFactory(
                    optd['submit_qtype'],
                    script,
                    cwd=optd['work_dir'],
                    run_time=43200,
                    name='benchmark',
                    max_array_size=optd['nproc'],
                    queue=optd['submit_queue'],
                    shell="/bin/bash",
            ) as task:
                task.run()
                task.wait(interval=5, monitor=monitor)

            # queue finished so unpickle results
            optd.update(ample_util.read_amoptd(optd['results_path']))
        else:
            benchmark_util.analyse(optd)
            ample_util.save_amoptd(optd)
        return

    @staticmethod
    def cleanup(optd):
        """Remove directories based on purge level
        """
        purge_level = optd['purge']
        to_remove = {
            1: ['modelling_workdir', 'ensembles_workdir'],
            2: ['models_dir', 'ensembles_directory', 'contacts_dir'],
        }
        for level in to_remove.keys():
            if purge_level >= level:
                for wdir in to_remove[level]:
                    if wdir in optd and optd[wdir] and os.path.isdir(optd[wdir]):
                        shutil.rmtree(optd[wdir])
        if purge_level >= 2:
            mrbump_util.purge_MRBUMP(optd)
        return

    def handle_contacts(self, optd):
        if optd["use_contacts"] and not optd['restraints_file']:
            con_util = contact_util.ContactUtil(
                optd['fasta'],
                'fasta',
                contact_file=optd['contact_file'],
                contact_format=optd['contact_format'],
                bbcontacts_file=optd['bbcontacts_file'],
                bbcontacts_format=optd["bbcontacts_format"],
                cutoff_factor=optd['restraints_factor'],
                distance_to_neighbor=optd['distance_to_neighbour'],
            )

            optd["contacts_dir"] = os.path.join(optd["work_dir"], "contacts")
            if not os.path.isdir(optd["contacts_dir"]):
                os.mkdir(optd["contacts_dir"])
            if con_util.require_contact_prediction:
                if con_util.found_ccmpred_contact_prediction_deps:
                    con_util.predict_contacts_from_sequence(wdir=optd["contacts_dir"])
                    optd["contact_file"] = con_util.contact_file
                    optd["contact_format"] = con_util.contact_format

            if con_util.do_contact_analysis:
                plot_file = os.path.join(optd['contacts_dir'], optd['name'] + ".cm.png")
                if optd['native_pdb'] and optd['native_pdb_std']:
                    structure_file = optd['native_pdb_std']
                elif optd["native_pdb"]:
                    structure_file = optd['native_std']
                else:
                    structure_file = None
                optd['contact_map'], optd['contact_ppv'] = con_util.summarize(
                    plot_file, structure_file, 'pdb', optd['native_cutoff']
                )
                restraints_file = os.path.join(optd['contacts_dir'], optd['name'] + ".cst")
                optd['restraints_file'] = con_util.write_restraints(
                    restraints_file, optd['restraints_format'], optd['energy_function']
                )
            else:
                con_util = None
        else:
            con_util = None
        return con_util

    def modelling_required(self, optd):
        return optd['make_frags'] or optd['make_models'] or optd['nmr_remodel']

    def ensembling(self, optd):
        if optd['import_ensembles']:
            ensembler.import_ensembles(optd)
        elif optd['ideal_helices']:
            ample_util.ideal_helices(optd)
            logger.info("*** Using ideal helices to solve structure ***")
            logger.warning('If ideal helices do not solve the structure, you may want to use -helical_ensembles in '
                           'place of -ideal_helices. AMPLE will then use a new set of helical ensembles which has been '
                           'very successful on solving challenging cases!')
        elif optd['helical_ensembles']:
            ample_util.ideal_helices(optd, ensembles=True)
            logger.info("*** Using helical ensembles to solve structure ***")
        else:
            # Check we have some models to work with
            if not (optd['single_model_mode'] or optd['processed_models']):
                ample_util.save_amoptd(optd)
                msg = "ERROR! Cannot find any pdb files in: {0}".format(optd['models_dir'])
                exit_util.exit_error(msg)
            optd['ensemble_ok'] = os.path.join(optd['work_dir'], 'ensemble.ok')
            if optd['submit_cluster']:
                # Pickle dictionary so it can be opened by the job to get the parameters
                ample_util.save_amoptd(optd)
                script = ensembler.cluster_script(optd)
                ensembler_timeout = ensembler.get_ensembler_timeout(optd)
                with TaskFactory(
                        optd['submit_qtype'],
                        script,
                        cwd=optd['work_dir'],
                        run_time=ensembler_timeout,
                        name='benchmark',
                        max_array_size=optd['nproc'],
                        queue=optd['submit_queue'],
                        shell="/bin/bash",
                ) as task:
                    task.run()
                    task.wait(interval=5, monitor=monitor)
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
                ensemble_summary = ensembler.ensemble_summary(optd['ensembles_data'])
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
        if not (optd['make_frags'] or optd['make_models'] or optd['nmr_remodel']):
            return
        # Set the direcotry where the final models will end up
        optd['models_dir'] = os.path.join(optd['work_dir'], 'models')
        if not os.path.isdir(optd['models_dir']):
            os.mkdir(optd['models_dir'])
        if not rosetta_modeller:
            rosetta_modeller = options_processor.process_rosetta_options(optd)
        # Make Rosetta fragments
        if optd['make_frags']:
            rosetta_modeller.generate_fragments(optd)
            optd['frags_3mers'] = rosetta_modeller.frags_3mers
            optd['frags_9mers'] = rosetta_modeller.frags_9mers
            optd['psipred_ss2'] = rosetta_modeller.psipred_ss2

        con_util = self.handle_contacts(optd)
        if optd['restraints_file']:
            rosetta_modeller.restraints_file = optd['restraints_file']

        if optd['make_models']:
            logger.info('----- making Rosetta models--------')
            logger.info('Making %s models...', optd['nmodels'])
            try:
                optd['processed_models'] = rosetta_modeller.ab_initio_model(processed_models=optd['processed_models'])
            except Exception as e:
                msg = "Error running ROSETTA to create models: {0}".format(e)
                exit_util.exit_error(msg, sys.exc_info()[2])
            logger.info('Modelling complete - models stored in: %s', optd['models_dir'])

        # Sub-select the decoys using contact information
        if con_util and optd['subselect_mode'] and not (optd['nmr_model_in'] or optd['nmr_remodel']):
            logger.info('Subselecting models from directory using provided contact information')
            subselect_data = con_util.subselect_decoys(
                optd['processed_models'], 'pdb', mode=optd['subselect_mode'], **optd
            )
            optd['processed_models'] = zip(*subselect_data)[0]
            optd['subselect_data'] = dict(subselect_data)

    def molecular_replacement(self, optd):
        mrbump_util.set_success_criteria(optd)
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
            logger.info("Running MRBUMP jobs in directory: %s", bump_dir)

            # Set an ensemble-specific phaser_rms if required
            if optd['phaser_rms'] == 'auto':
                ensembler.set_phaser_rms_from_subcluster_score(optd)

            # Sort the ensembles in a favourable way
            logger.info("Sorting ensembles")
            sort_keys = ['cluster_num', 'truncation_level', 'subcluster_radius_threshold', 'side_chain_treatment']
            ensemble_pdbs_sorted = ensembler.sort_ensembles(
                optd['ensembles'], optd['ensembles_data'], keys=sort_keys, prioritise=True
            )

            # Create job scripts
            logger.info("Generating MRBUMP runscripts")
            optd['mrbump_scripts'] = mrbump_util.write_mrbump_files(
                ensemble_pdbs_sorted,
                optd,
                job_time=mrbump_util.MRBUMP_RUNTIME,
                ensemble_options=optd['ensemble_options'],
                directory=bump_dir,
            )

        # Create function for monitoring jobs - static function decorator?
        if self.ample_output:

            def monitor():
                r = mrbump_util.ResultsSummary()
                r.extractResults(optd['mrbump_dir'], purge=bool(optd['purge']))
                optd['mrbump_results'] = r.results
                return self.ample_output.display_results(optd)

        else:
            monitor = None

        # Save results here so that we have the list of scripts and mrbump directory set
        ample_util.save_amoptd(optd)

        # Change to mrbump directory before running
        os.chdir(optd['mrbump_dir'])

        processes = optd['nproc']
        submit_max_array = optd['submit_max_array']
        if optd['submit_qtype'] != 'local':
            processes = 1
            submit_max_array = optd['nproc']

        with TaskFactory(
                optd['submit_qtype'],
                optd['mrbump_scripts'],
                cwd=bump_dir,
                run_time=mrbump_util.MRBUMP_RUNTIME,
                name="mrbump",
                processes=processes,
                max_array_size=submit_max_array,
                queue=optd['submit_queue'],
                shell="/bin/bash",
        ) as task:
            task.run()

            if optd['early_terminate']:
                task.wait(interval=5, monitor=monitor, success_f=mrbump_util.checkSuccess)
            else:
                task.wait(interval=5, monitor=monitor)

        if not task.completed:
            msg = (
                    "An error code was returned after running MRBUMP on the ensembles!\n"
                    + "For further information check the logs in directory: {0}".format(optd['mrbump_dir'])
            )
            logger.critical(msg)

        # Collect the MRBUMP results
        results_summary = mrbump_util.ResultsSummary()
        optd['mrbump_results'] = results_summary.extractResults(optd['mrbump_dir'], purge=bool(optd['purge']))
        optd['success'] = results_summary.success
        ample_util.save_amoptd(optd)
        summary = mrbump_util.finalSummary(optd)
        logger.info(summary)

    def process_models(self, optd):
        process_models.extract_and_validate_models(optd)
        # Need to check if Quark and handle things accordingly
        if optd['quark_models']:
            # We always add sidechains to QUARK models if SCWRL is installed
            if ample_util.is_exe(optd['scwrl_exe']):
                optd['use_scwrl'] = True
            else:
                # No SCWRL so don't do owt with the side chains
                logger.info('Using QUARK models but SCWRL is not installed ' 'so only using %s sidechains', UNMODIFIED)
                optd['side_chain_treatments'] = [UNMODIFIED]
        ample_util.save_amoptd(optd)

    def setup(self, optd):
        """We take and return an ample dictionary as an argument.

        This is required because options_processor.process_restart_options Changes what
        optd points at, and this means that if we just use the reference, we end up
        pointing at the old, obsolete dictionary

        """
        optd = options_processor.restart_amoptd(optd)
        optd['ccp4_version'] = ample_util.CCP4.version.version
        logger.info(reference_manager.header)
        logger.info("AMPLE version: %s", str(version.__version__))
        logger.info("Using CCP4 version: %s from directory: %s", ample_util.CCP4.version, ample_util.CCP4.root)
        logger.info("Running on host: %s", platform.node())
        logger.info("Running on platform: %s", platform.platform())
        logger.info('Running on %d processors', optd['nproc'])
        logger.info("Job started at: %s", time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
        logger.info("Invoked with command-line:\n%s\n", " ".join(sys.argv))
        logger.info("Running in directory: %s\n", optd['work_dir'])
        if pyrvapi_results.pyrvapi:
            self.ample_output = pyrvapi_results.AmpleOutput(optd)
            self.ample_output.display_results(optd)
        options_processor.check_mandatory_options(optd)
        optd = options_processor.process_restart_options(optd)
        if not optd['restart_pkl']:
            options_processor.process_options(optd)
        if optd['dry_run']:
            logger.info('Dry run finished checking options - cleaning up...')
            os.chdir(optd['run_dir'])
            shutil.rmtree(optd['work_dir'])
            sys.exit(0)
        logger.info('All needed programs are found, continuing...')
        return optd

    def setup_workdir(self, argso):
        """Make a work directory - this way all output goes into this directory.
        
        This is done before the loggers has been set up so no logging is possible.
        """
        if argso['work_dir'] and not argso['restart_pkl']:
            try:
                os.mkdir(argso['work_dir'])
            except Exception as e:
                msg = "Cannot create work_dir {0}: {1}".format(argso['work_dir'], e)
                exit_util.exit_error(msg, sys.exc_info()[2])
        if not argso['work_dir']:
            if not os.path.exists(argso['run_dir']):
                msg = 'Cannot find run directory: {0}'.format(argso['run_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])
            if argso['rvapi_document']:
                # With JSCOFE we run in the run directory
                argso['work_dir'] = argso['run_dir']
            else:
                argso['work_dir'] = ample_util.make_workdir(argso['run_dir'], ccp4i2=bool(argso['ccp4i2_xml']))
        os.chdir(argso['work_dir'])
        return argso['work_dir']


if __name__ == "__main__":
    try:
        Ample().main()
    except Exception as e:
        msg = "Error running main AMPLE program: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
