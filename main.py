#!/usr/bin/env ccp4-python
"""
This is AMPLE

This script is named ample.py due to a problem with running the multiprocessing 
(which is used to parallelise the running of jobs on a local machine - see python/workers.py) 
module under windows. The multiprocessing module on windows requires that it 
can import the main module, and the import machinery requires that any file 
being imported is named <foo>.py, and any changes to this would require hacking 
the multiprocessing module, so to avoid this, our script must be called ample.py
"""

# python imports
import cPickle
import glob
import logging
import os
import platform
import shutil
import sys
import time

# Our imports
from ample.modelling import rosetta_model
from ample.python import pyrvapi_results
from ample.util import ample_util
from ample.util import argparse_util
from ample.util import benchmark_util
from ample.util import config_util
from ample.util import contacts_util
from ample.util import ensembler_util
from ample.util import exit_util
from ample.util import mrbump_util
from ample.util import mtz_util
from ample.util import options_util
from ample.util import options_processor
from ample.util import pdb_edit
from ample.util import sequence_util
from ample.util import workers_util
from ample.util import version


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
    formatter = logging.Formatter('%(asctime)s - %(name)s [%(lineno)d] - %(levelname)s - %(message)s')
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
            msg += "Cannot find the $CCP4_SCR directory: {0}\n".format(os.environ['CCP4_SCR'])
            msg += "The directory will be created, but it should have already been created by the CCP4 startup scripts\n"
            msg += "Please make sure CCP4 is installed and the setup scripts have been run."
            logger.critical(msg)
            os.mkdir(os.environ['CCP4_SCR'])
            #exit_util.exit_error(msg)
    
        # Record the CCP4 version we're running with  - also required in pyrvapi_results
        amoptd['ccp4_version'] = ample_util.ccp4_version()
        
        return os.environ['CCP4']

    def run(self, args=None):
        """Main AMPLE routine.
        
        We require this as the multiprocessing module (only on **!!*%$$!! Windoze) 
        requires that the main module can be imported. We there need ample to be 
        a python script that can be imported, hence the main routine with its 
        calling protected by the if __name__=="__main__":...
        
        args is an option argument that can contain the command-line arguments 
        for the program - required for testing.
        """ 
        argso = argparse_util.process_command_line(args=args)
        self.amopt = amopt = options_util.AmpleOptions()
        amopt.populate(argso)
        
        # Make a work directory - this way all output goes into this directory
        if amopt.d['work_dir']:
            logger.info('Making a named work directory: {0}'.format(amopt.d['work_dir']))
            try:
                os.mkdir(amopt.d['work_dir'])
            except:
                msg = "Cannot create work_dir {0}".format(amopt.d['work_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])
        else:
            if not os.path.exists(amopt.d['run_dir']):
                msg = 'Cannot find run directory: {0}'.format(amopt.d['run_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])
            logger.info('Making a run directory: checking for previous runs...')
            amopt.d['work_dir'] = ample_util.make_workdir(amopt.d['run_dir'], 
                                                          ccp4_jobid=amopt.d['ccp4_jobid'])
        # Go to the work directory
        os.chdir(amopt.d['work_dir'])
        
        # Set up logging
        ample_log = os.path.join(amopt.d['work_dir'], 'AMPLE.log')
        debug_log = os.path.join(amopt.d['work_dir'], 'debug.log')
        amopt.d['ample_log'] = ample_log
        
        setup_file_logging(ample_log, debug_log)
        
        # Make sure the CCP4 environment is set up properly
        ccp4_home = self.setup_ccp4(amopt.d)
        ccp4_version = ".".join([str(x) for x in amopt.d['ccp4_version']])
        
        # Print out Version and invocation
        logger.info(ample_util.header)
        logger.info("AMPLE version: {0}".format(version.__version__))
        logger.info("Running with CCP4 version: {0} from directory: {1}".format(ccp4_version, ccp4_home))
        logger.info("Job started at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        logger.info("Running on host: {0}".format(platform.node()))
        logger.info("Invoked with command-line:\n{0}\n".format(" ".join(sys.argv)))
        logger.info("Running in directory: {0}\n".format(amopt.d['work_dir']))
        
        # Display pyrvapi results
        if pyrvapi_results.pyrvapi:
            self.output_gui = pyrvapi_results.AmpleOutput()
            self.output_gui.display_results(amopt.d)
        
        # Check mandatory/exclusive options
        options_processor.check_mandatory_options(amopt.d)
        
        # Check if we are restarting from an existing pkl file - we don't process the options from this
        # run if so
        amopt.d = options_processor.process_restart_options(amopt.d)
        if not amopt.d['restart_pkl']:
            options_processor.process_options(amopt.d) # Only process the remaining options if we aren't in restart mode
        rosetta_modeller = options_processor.process_rosetta_options(amopt.d)
        
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
            
        if amopt.d['benchmark_mode'] and amopt.d['native_pdb']:
            # Process the native before we do anything else
            benchmark_util.analysePdb(amopt.d)       
    
        # Make Rosetta fragments
        if amopt.d['make_frags']:
            rosetta_modeller.generate_fragments(amopt.d)
            amopt.d['frags_3mers'] = rosetta_modeller.frags_3mers
            amopt.d['frags_9mers'] = rosetta_modeller.frags_9mers
            amopt.d['psipred_ss2'] = rosetta_modeller.psipred_ss2
    
        # In case file created above we need to tell the rosetta_modeller where it is
        # otherwise not used as not created before object initialised    
        if amopt.d['make_models'] and (amopt.d['use_contacts'] or amopt.d['restraints_file']):
            cm = contacts_util.Contacter(optd=amopt.d)
            
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
                    exit_util.exit_error(msg, sys.exc_info()[2])
            else:
                logger.info('making {0} models...'.format(amopt.d['nmodels']))
                try:
                    rosetta_modeller.ab_initio_model(monitor=monitor)
                except Exception, e:
                    msg = "Error running ROSETTA to create models: {0}".format(e)
                    exit_util.exit_error(msg, sys.exc_info()[2])
                if not pdb_edit.check_pdb_directory(amopt.d['models_dir'], sequence=amopt.d['sequence']):
                    msg = "Problem with rosetta pdb files - please check the log for more information"
                    exit_util.exit_error(msg)
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
                        logger.info('Using QUARK models but SCWRL is not installed so only using {0} sidechains'.format(UNMODIFIED))
                        amopt.d['side_chain_treatments'] = [ UNMODIFIED ]
    
        # Save the results
        ample_util.saveAmoptd(amopt.d)
        
        if amopt.d['make_ensembles']:
            if amopt.d['import_ensembles']:
                ensembler_util.import_ensembles(amopt.d)
            elif amopt.d['ideal_helices']:
                amopt.d['ensembles'], amopt.d['ensemble_options'], amopt.d['ensembles_data'] = ample_util.ideal_helices(amopt.d['fasta_length'])
                logger.info("*** Using ideal helices to solve structure ***")
            else:
                # Check we have some models to work with
                if not (amopt.d['cluster_method'] is 'import' or amopt.d['single_model_mode']) and \
                   not glob.glob(os.path.join(amopt.d['models_dir'], "*.pdb")):
                    ample_util.saveAmoptd(amopt.d)
                    msg = "ERROR! Cannot find any pdb files in: {0}".format(amopt.d['models_dir'])
                    exit_util.exit_error(msg)
                amopt.d['ensemble_ok'] = os.path.join(amopt.d['work_dir'],'ensemble.ok')
                if amopt.d['submit_cluster']:
                    # Pickle dictionary so it can be opened by the job to get the parameters
                    ample_util.saveAmoptd(amopt.d)
                    script = ensembler_util.cluster_script(amopt.d)
                    workers_util.run_scripts(job_scripts=[script],
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
                    try: ensembler_util.create_ensembles(amopt.d)
                    except Exception, e:
                        msg = "Error creating ensembles: {0}".format(e)
                        exit_util.exit_error(msg, sys.exc_info()[2])
                        
                # Check we have something to work with
                if not os.path.isfile(amopt.d['ensemble_ok']) or not amopt.d.has_key('ensembles') or not len(amopt.d['ensembles']):
                    msg = "Problem generating ensembles!"
                    exit_util.exit_error(msg)
                    
                if not (amopt.d['homologs'] or amopt.d['single_model_mode']):
                    ensemble_summary = ensembler_util.ensemble_summary(amopt.d['ensembles_data'])
                    logger.info(ensemble_summary)
                
            # Save the results
            ample_util.saveAmoptd(amopt.d)
            
            # Bail here if we didn't create anything
            if not len(amopt.d['ensembles']):
                msg = "### AMPLE FAILED TO GENERATE ANY ENSEMBLES! ###\nExiting..."
                exit_util.exit_error(msg)
        
        # Update results view
        if self.output_gui: self.output_gui.display_results(amopt.d)
         
        if amopt.d['do_mr']:
            if not amopt.d['mrbump_scripts']:
                # MRBUMP analysis of the ensembles
                logger.info('----- Running MRBUMP on ensembles--------\n\n')
                if len(amopt.d['ensembles']) < 1:
                    msg = "ERROR! Cannot run MRBUMP as there are no ensembles!"
                    exit_util.exit_error(msg)
                 
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
                ensemble_pdbs_sorted = ensembler_util.sort_ensembles(amopt.d['ensembles'],
                                                               amopt.d['ensembles_data'])

                # Create job scripts
                logger.info("Generating MRBUMP runscripts")
                amopt.d['mrbump_scripts'] = mrbump_util.write_mrbump_files(ensemble_pdbs_sorted,
                                                                            amopt.d,
                                                                            job_time=mrbump_util.MRBUMP_RUNTIME,
                                                                            ensemble_options=amopt.d['ensemble_options'],
                                                                            directory=bump_dir )
            # Create function for monitoring jobs - static function decorator?
            if self.output_gui:
                def monitor():
                    r = mrbump_util.ResultsSummary()
                    r.extractResults(amopt.d['mrbump_dir'], purge=amopt.d['purge'])
                    amopt.d['mrbump_results'] = r.results
                    return self.output_gui.display_results(amopt.d)
            else:
                monitor = None
                
            # Save results here so that we have the list of scripts and mrbump directory set
            ample_util.saveAmoptd(amopt.d)
                
            # Change to mrbump directory before running
            os.chdir(amopt.d['mrbump_dir'])  
            ok = workers_util.run_scripts(job_scripts=amopt.d['mrbump_scripts'],
                                          monitor=monitor,
                                          check_success=mrbump_util.checkSuccess,
                                          early_terminate=amopt.d['early_terminate'],
                                          chdir=False,
                                          nproc=amopt.d['nproc'],
                                          job_time=mrbump_util.MRBUMP_RUNTIME,
                                          job_name='mrbump',
                                          submit_cluster=amopt.d['submit_cluster'],
                                          submit_qtype=amopt.d['submit_qtype'],
                                          submit_queue=amopt.d['submit_queue'],
                                          submit_array=amopt.d['submit_array'],
                                          submit_max_array=amopt.d['submit_max_array'])
         
            if not ok:
                msg = "Error running MRBUMP on the ensembles!\nCheck logs in directory: {0}".format(amopt.d['mrbump_dir'])
                exit_util.exit_error(msg)
        
            # Collect the MRBUMP results
            results_summary = mrbump_util.ResultsSummary()
            amopt.d['mrbump_results'] = results_summary.extractResults(amopt.d['mrbump_dir'], purge=amopt.d['purge'])
            amopt.d['success'] = results_summary.success
            
            ample_util.saveAmoptd(amopt.d)
        
            # Now print out the final summary
            summary = mrbump_util.finalSummary(amopt.d)
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
                script = benchmark_util.cluster_script(amopt.d)
                workers_util.run_scripts(job_scripts=[script],
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
                benchmark_util.analyse(amopt.d)
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
        exit_util.exit_error(msg, sys.exc_info()[2])
