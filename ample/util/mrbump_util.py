#!/usr/bin/env ccp4-python

import copy
import pickle
import glob
import logging
import os
import shutil
import sys

# Hack to make sure we can find the modules we need
if __name__ == "__main__":
    root = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2])
    sys.path.insert(0, os.path.join(root, "scripts"))

from ample.util import ample_util
from ample.util import mrbump_cmd
from ample.util import printTable

# FS [15/11/2018] -> New MRBUMP structure allows direct imports, wait for release then replace
mrbumpd = os.path.join(os.environ['CCP4'], "share", "mrbump", "include", "parsers")
sys.path.insert(0, mrbumpd)
import parse_arpwarp
import parse_buccaneer
import parse_phaser

TOP_KEEP = 3 # How many of the top shelxe/phaser results to keep for the gui
MRBUMP_RUNTIME = 172800 # allow 48 hours for each mrbump job
SHELXE_MAX_PERMITTED_RESOLUTION = 3.0

# Values to determine when job has succeeded - required at module level this may be set by AMPLE from
# the command line
SUCCESS_PHASER_TFZ = 8.0
SUCCESS_PHASER_LLG = 120
SUCCESS_RFREE = 0.4
SUCCESS_SHELXE_CC = 25.0
SUCCESS_SHELXE_ACL = 10


# We need a null logger so that we can be used without requiring a logger
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())


class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """

    def __init__(self, results=None, results_pkl=None):
        """
        Parameters
        ----------
        results_pkl : file
           A pickled AMPLE results dictionary
        """ 
        self.results = []
        self.pname = "archive"
        self.pdir = None
        self.success = False
        if results_pkl and os.path.isfile(results_pkl):
            with open(results_pkl) as f:
                resd = pickle.load(f)
            mkey = 'mrbump_results'
            if mkey in resd and len(resd[mkey]):
                self.results = resd[mkey]
        elif results:
            self.results = results
        return

    def analyseResult(self, result):
        mrDir = result["MR_directory"]
        # result.ensembleName = result.name[9:-6]
        if result["MR_program"] == "PHASER":
            if result["PHASER_pdbout"]:
                phaserP = parse_phaser.PhaserPdbParser(result["PHASER_pdbout"])
                result["PHASER_LLG"] = phaserP.LLG
                result["PHASER_TFZ"] = phaserP.TFZ
            
            phaserLog = os.path.join(mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(result["MR_program"].lower(), result['ensemble_name']))
            if os.path.isfile(phaserLog):
                phaserP = parse_phaser.PhaserLogParser(phaserLog, noLLG=True)
                # result.phaserLog    = phaserLog
                result["PHASER_time"] = phaserP.time
                result["PHASER_killed"] = phaserP.killed
        buccaneerLog = os.path.join(mrDir,
                                     "build/shelxe/rebuild/buccaneer",
                                     "buccaneer.log")
        bp = parse_buccaneer.BuccaneerLogParser()
        if os.path.isfile(buccaneerLog):
            bp.parse(buccaneerLog)
            result["SXRBUCC_final_Rfree"] = bp.finalRfree
            result["SXRBUCC_final_Rfact"] = bp.finalRfact
        arpLog = os.path.join(mrDir,
                              "build/shelxe/rebuild/arpwarp",
                              "arpwarp.log")
        if os.path.isfile(arpLog):
            ap = parse_arpwarp.ArpwarpLogParser()
            ap.parse(arpLog)
            result["SXRARP_final_Rfact"] = ap.finalRfact
            result["SXRARP_final_Rfree"] = ap.finalRfree
        return



    def createDict(self):
        d = {}
        # our additional keys
        d['ensemble_name'] = None
        d['MR_program'] = None
        d['name'] = None
        d['Search_directory'] = None
        d['MR_directory'] = None
        d['Solution_Type'] = None
        
        d['PHASER_LLG'] = None
        d['PHASER_TFZ'] = None
        d['PHASER_RFZ'] = None
        d['PHASER_time'] = None
        d['PHASER_killed'] = None
        d['PHASER_pdbout'] = None
        d['PHASER_mtzout'] = None
        d['PHASER_logfile'] = None
        d['PHASER_version'] = None
        d['PHASER_error'] = None
        
        d['MOLREP_score'] = None
        d['MOLREP_time'] = None
        d['MOLREP_pdbout'] = None
        d['MOLREP_logfile'] = None
        d['MOLREP_version'] = None
        
        d['REFMAC_Rfact'] = None
        d['REFMAC_Rfree'] = None
        d['REFMAC_pdbout'] = None
        d['REFMAC_mtzout'] = None
        d['REFMAC_logfile'] = None
        d['REFMAC_version'] = None
        
        d['BUCC_final_Rfact'] = None
        d['BUCC_final_Rfree'] = None
        d['BUCC_pdbout'] = None
        d['BUCC_mtzout'] = None
        d['BUCC_logfile'] = None
        d['BUCC_version'] = None
        
        d['ARP_final_Rfact'] = None
        d['ARP_final_Rfree'] = None
        d['ARP_pdbout'] = None
        d['ARP_mtzout'] = None
        d['ARP_logfile'] = None
        d['ARP_version'] = None
        
        d['SHELXE_CC'] = None
        d['SHELXE_ACL'] = None
        d['SHELXE_MCL'] = None
        d['SHELXE_NC'] = None
        d['SHELXE_wMPE'] = None
        d['SHELXE_os'] = None
        d['SHELXE_time'] = None
        d['SHELXE_pdbout'] = None
        d['SHELXE_phsout'] = None
        d['SHELXE_mtzout'] = None
        d['SHELXE_logfile'] = None
        d['SHELXE_version'] = None
        
        d['SXRBUCC_version'] = None
        d['SXRBUCC_final_Rfact'] = None
        d['SXRBUCC_final_Rfree'] = None
        d['SXRBUCC_pdbout'] = None
        d['SXRBUCC_mtzout'] = None
        d['SXRBUCC_logfile'] = None
        
        d['SXRARP_version'] = None
        d['SXRARP_final_Rfact'] = None
        d['SXRARP_final_Rfree'] = None
        d['SXRARP_pdbout'] = None
        d['SXRARP_mtzout'] = None
        d['SXRARP_logfile'] = None
        return d
    
    def _extractPurged(self, mrbump_dir):
        """Recreate a list of the jobs that have been purged"""
        purged_results = {}
        self.pdir = os.path.join(mrbump_dir, self.pname)
        if not os.path.isdir(self.pdir):
            os.mkdir(self.pdir)
        pkls = glob.glob(os.path.join(self.pdir, "*.pkl"))
        if pkls:
            for p in pkls:
                with open(p) as f:
                    d = pickle.load(f)
                purged_results[d['ensemble_name']] = d
        return purged_results

    def extractResults(self, mrbump_dir, purge=False):
        if not mrbump_dir or not os.path.isdir(mrbump_dir):
            raise RuntimeError("Cannot find mrbump_dir: {0}".format(mrbump_dir))
        purged_results = {}
        if purge:
            purged_results = self._extractPurged(mrbump_dir)
        self._extractResults(mrbump_dir, archived_ensembles=purged_results.keys())
        if purge:
            self._purgeFailed()
            self.results += purged_results.values()
        self.sortResults()
        self.success = any([self.jobSucceeded(r) for r in self.results])
        return self.results

    def _extractResults(self, mrbump_dir, archived_ensembles=None):
        """
        Find the results from running MRBUMP and sort them
        """
        mrbump_dir = os.path.abspath(mrbump_dir)
        if not os.path.isdir(mrbump_dir):
            logger.warn("extractResults - is not a valid directory: {0}".format(mrbump_dir))
            return []
        # Get a list of the ensembles (could get this from the amopt dictionary)
        # For now we just use the submission scripts and assume all have .sh or .sub extension
        ext = '.sh'
        if sys.platform.startswith("win"):
            ext = '.bat'
        ensembles = [ os.path.splitext(os.path.basename(e))[0] for e in glob.glob(os.path.join(mrbump_dir, "*" + ext))]
        if not len(ensembles):
            # legacy - try .sub
            ensembles = [ os.path.splitext(os.path.basename(e))[0] for e in glob.glob(os.path.join(mrbump_dir, "*.sub"))]
        if not len(ensembles):
            logger.warn("Could not extract any results from directory: {0}".format(mrbump_dir))
            return []
        # reset any results
        results = []
        failed = {}  # dict mapping failures to what went wrong - need to process at the end
        for ensemble in ensembles:
            # Skip ones that we've archived
            if archived_ensembles and ensemble in archived_ensembles:
                continue
            # Check job directory
            jobDir = os.path.join(mrbump_dir, 'search_' + ensemble + '_mrbump')
            if not os.path.isdir(jobDir):
                jobDir = os.path.join( mrbump_dir, 'search_'+ensemble )
            if not os.path.isdir(jobDir):
                # As we call this every time we monitor a job running, we don't want to print this out all the time
                # logger.debug("Missing job directory: {0}".format(jobDir))
                failed[ensemble] = "no_job_directory"
                continue
            logger.debug(" -- checking directory for results: {0}".format(jobDir))
            # Check if finished
            if not os.path.exists(os.path.join(jobDir, "results", "finished.txt")):
                logger.debug("Found unfinished job: {0}".format(jobDir))
                failed[ensemble] = "unfinished"
                continue
            # Check resultsTable.dat
            resultsDict = os.path.join(jobDir, "results", "resultsTable.pkl")
            if os.path.isfile(resultsDict):
                results += self.processMrbumpPkl(resultsDict)
            else:
                logger.debug(" -- Could not find results files: {0}".format(resultsDict))
                failed[ensemble] = "missing-results-file"
                continue
        # Process the failed results
        if failed:
            results += self._processFailed(mrbump_dir, failed)
        if not len(results):
            logger.warn("Could not extract any results from directory: {0}".format(mrbump_dir))
        self.results = results
        return

    @staticmethod
    def jobSucceeded(job_dict):
        success = False
        if 'SHELXE_CC' in job_dict and job_dict['SHELXE_CC'] and float(job_dict['SHELXE_CC']) >= SUCCESS_SHELXE_CC and \
           'SHELXE_ACL' in job_dict and job_dict['SHELXE_ACL'] and float(job_dict['SHELXE_ACL']) >= SUCCESS_SHELXE_ACL:
            success = True
        elif 'BUCC_final_Rfree' in job_dict and job_dict['BUCC_final_Rfree'] and float(job_dict['BUCC_final_Rfree']) <= SUCCESS_RFREE:
            success = True
        elif 'ARP_final_Rfree' in job_dict and job_dict['ARP_final_Rfree'] and float(job_dict['ARP_final_Rfree']) <= SUCCESS_RFREE:
            success = True
        elif 'REFMAC_Rfree' in job_dict and job_dict['REFMAC_Rfree'] and float(job_dict['REFMAC_Rfree']) <= SUCCESS_RFREE:
            success = True
        elif 'PHASER_LLG' in job_dict and 'PHASER_TFZ' in job_dict and job_dict['PHASER_LLG'] and job_dict['PHASER_TFZ'] and \
        float(job_dict['PHASER_LLG']) >= SUCCESS_PHASER_LLG and float(job_dict['PHASER_TFZ']) >= SUCCESS_PHASER_TFZ:
            success = True
        return success

    def processMrbumpPkl(self, resultsPkl):
        """Process dictionary
        """
        with open(resultsPkl) as f:
            rD = pickle.load(f)
        if not rD:
            return []
        results = []
        for name, d1 in rD.iteritems():
            for mrprog, d2 in d1.iteritems():
                # Check if all the entries are None - means job didn't run.
                # Should probably think of a better way to spot that (Search_directory is always set)
                if not any([v for k, v in d2.iteritems() if k != 'Search_directory']):
                    continue
                # Add MR program as dictionary entry
                d = copy.copy(d2)
                del d['SearchModel_filename']
                d['name'] = name
                # name is e.g.: loc0_ALL_c1_tl100_r2_allatom_UNMOD
                d['ensemble_name'] = name[9:-6]
                d['MR_program'] = mrprog
                # Hack for old versions
                if 'JobDirectory' in d:
                    d['MR_directory'] = d['JobDirectory']
                    del d['JobDirectory']
                    d['Search_directory'] = os.sep.join(d['MR_directory'].split(os.sep)[:-5])
                if 'final_Rfree' in d:
                    d['REFMAC_Rfree'] = d['final_Rfree']
                    d['REFMAC_Rfact'] = d['final_Rfact']
                    del d['final_Rfree']
                    del d['final_Rfact']
                results.append(d)
        return results
 
    def _processFailed(self, mrbump_dir, failed):
        """Generate dictionaries for failed results
        """
        results = []
        for ensemble, reason in failed.iteritems():
            d = self.createDict()
            # name hard-coded
            # d['name'] = "loc0_ALL_" + ensemble + "_UNMOD"
            d['name'] = "loc0_ALL_" + ensemble + "_UNMOD"
            d['ensemble_name'] = ensemble 
            d['Search_directory'] = os.path.join(mrbump_dir, 'search_' + ensemble + '_mrbump')
            d['Solution_Type'] = reason
            results.append(d)
        logger.debug("Added {0} MRBUMP result failures".format(len(failed)))
        return results
    
    def _purgeFailed(self):
        """Remove the MRBUMP directories of any jobs that don't pass the keep criteria and archive their job dictionaries"""
        completed = [r for r in self.results if not(job_unfinished(r))]
        if completed:
            to_keep = []
            min_len = min(len(completed), TOP_KEEP)
            for r in ResultsSummary.sortResultsStatic(completed, prioritise='SHELXE_CC')[:min_len]:
                if r not in to_keep:
                    to_keep.append(r)
            for r in ResultsSummary.sortResultsStatic(completed, prioritise='PHASER_TFZ')[:min_len]:
                if r not in to_keep:
                    to_keep.append(r)
            for r in completed:
                if r not in to_keep:
                    pkl = os.path.join(self.pdir, "{0}.pkl".format(r['ensemble_name']))
                    with open(pkl, 'w') as f:
                        pickle.dump(r, f)
                    shutil.rmtree(r['Search_directory'])
        
    def results_table(self, results):
        resultsTable = []
        keys = ['ensemble_name', 'MR_program', 'Solution_Type']
        keys += _resultsKeys(results)
        resultsTable.append(keys)
        for r in results:
            resultsTable.append([r[k] for k in keys])
        return resultsTable
    

    def sortResults(self, prioritise=False):
        """Wrapper function to allow calls with self"""
        self.results = ResultsSummary.sortResultsStatic(self.results)

    @staticmethod
    def sortResultsStatic(results, prioritise=False):
        """Sort the results"""
        SHELXE = False
        BUCC = False
        ARP = False
        REFMAC = False
        PHASER = False
        for r in results:
            if 'SHELXE_CC' in r and r['SHELXE_CC'] and float(r['SHELXE_CC']) > 0.0:
                SHELXE = True
            if 'BUCC_final_Rfact' in r and r['BUCC_final_Rfact'] and float(r['BUCC_final_Rfact']) < 1.0:
                BUCC = True
            if 'ARP_final_Rfree' in r and r['ARP_final_Rfree'] and float(r['ARP_final_Rfree']) < 1.0:
                ARP = True
            if 'REFMAC_Rfree' in r and r['REFMAC_Rfree'] and float(r['REFMAC_Rfree']) < 1.0:
                REFMAC = True
            if 'PHASER_TFZ' in r and r['PHASER_TFZ'] and float(r['PHASER_TFZ']) > 0.0:
                PHASER = True
            
        reverse = False
        sortf = False
        if SHELXE and prioritise != "PHASER_TFZ":
            reverse = True
            sortf = lambda x: float(0) if x['SHELXE_CC']  is None else float(x['SHELXE_CC'])
        elif BUCC and not prioritise == "PHASER_TFZ":
            sortf = lambda x: float('inf') if x['BUCC_final_Rfact']  is None else float(x['BUCC_final_Rfact'])
        elif ARP and not prioritise == "PHASER_TFZ":
            sortf = lambda x: float('inf') if x['ARP_final_Rfree']  is None else float(x['ARP_final_Rfree'])
        elif REFMAC and not prioritise == "PHASER_TFZ":
            sortf = lambda x: float('inf') if x['REFMAC_Rfree']  is None else float(x['REFMAC_Rfree'])
        elif PHASER:
            reverse = True
            sortf = lambda x: float(0) if x['PHASER_TFZ']  is None else float(x['PHASER_TFZ'])
            
        if sortf:
            results.sort(key=sortf, reverse=reverse)
        return results

    def summariseResults(self, mrbump_dir):
        """Return a string summarising the results"""
        results = self.extractResults(mrbump_dir)
        if len(results):
            return self.summaryString()
        else:
            return "\n!!! No results found in directory: {0}\n".format(mrbump_dir)

    def summaryString(self):
        """Return a string suitable for printing the sorted results"""
        resultsTable = self.results_table(self.results)
        # Format the results
        table = printTable.Table()
        summary = table.pprint_table(resultsTable)
        r = "\n\nOverall Summary:\n\n"
        r += summary
        # Hack need to think of a better way to do this when there are no valid results
        top = self.results[0]
        k = None
        for p in ['Search_directory','MR_directory']:
            if p in top.keys():
                k = p
        assert k,"Missing search directory key in results dictionary"
        if top[k]:
            r += '\nBest Molecular Replacement results so far are in:\n\n'
            r += top[k]
        r += '\n\n'
        return r
    
    def topFiles(self, num_results=3):
        """Return a list of dictionaries listing the top num_results PDB and MTZ files
        
        Parameters
        ----------
        num_results : int
           How many of the top results to return
    
        Returns
        -------
        topf : list
           A list of dictionaries, one per result, with xyz, mtz and info keys
        
        """
        topf = []
        # list of PDB, MTZ, Explanation of file type - ordered by their desirability
        poss = [ ('SXRARP', 'SXRARP_pdbout','SXRARP_mtzout', 'ARPWARP rebuild of SHELXE trace of MR result'),
                 ('SXRBUCC', 'SXRBUCC_pdbout','SXRBUCC_mtzout', 'BUCCANEER rebuild of SHELXE trace of MR result'),
                 ('SHELXE', 'SHELXE_pdbout','SHELXE_mtzout', 'SHELXE trace of MR result'),
                 ('ARP', 'ARP_pdbout','ARP_mtzout', 'ARPWARP rebuild of MR result'),
                 ('BUCC', 'BUCC_pdbout','BUCC_mtzout', 'BUCCANEER rebuild of MR result'),
                 ('REFMAC,', 'REFMAC_pdbout','REFMAC_mtzout', 'REFMAC-refined MR result') ]
        for result in self.results[0 : min(num_results, len(self.results)+1) ]:
            for stype, pdb_key, mtz_key, source in poss:
                if pdb_key in result and result[pdb_key] and mtz_key in result and result[mtz_key]:
                    # Don't check paths for now as it screws up unittests as files don't actually exist
                    #if not (os.path.isfile(result[pdb_key]) and os.path.isfile(result[mtz_key])): continue
                    topf.append({ 'name' : result['ensemble_name'], 
                                  'type' : stype,
                                  'info' : source,
                                  'pdb' : result[pdb_key],
                                  'mtz' : result[mtz_key] })
                    break # Stop as soon as we find one
        if len(topf):
            return topf
#
# Module functions
#
def _resultsKeys(results):
    keys = []
    # Build up list of keys we want to print based on what we find in the results
    if any([True for r in results if r['PHASER_LLG']]):
        keys += ['PHASER_LLG']
    if any([True for r in results if r['PHASER_TFZ']]):
        keys += ['PHASER_TFZ']
    if any([True for r in results if r['REFMAC_Rfree'] and r['REFMAC_Rfree'] < 1.0]):
        keys += ['REFMAC_Rfact', 'REFMAC_Rfree']
    if any([True for r in results if r['BUCC_final_Rfact'] and r['BUCC_final_Rfact'] < 1.0]):
        keys += ['BUCC_final_Rfact', 'BUCC_final_Rfree']
    if any([True for r in results if r['ARP_final_Rfact'] and r['ARP_final_Rfact'] < 1.0]):
        keys += ['ARP_final_Rfact', 'ARP_final_Rfree']
    if any([True for r in results if r['SHELXE_CC']]):
        keys += ['SHELXE_CC', 'SHELXE_ACL']
    if any([True for r in results if r['SXRBUCC_final_Rfact']]):
        keys += ['SXRBUCC_final_Rfact', 'SXRBUCC_final_Rfree']
    if any([True for r in results if r['SXRARP_final_Rfact']]):
        keys += ['SXRARP_final_Rfact', 'SXRARP_final_Rfree']
    return keys


def checkSuccess(script_path):
    """Check if a job ran successfully.
    
    Parameters
    ----------
    script_path : str
       Path to the MrBUMP script
    
    Returns
    -------
    bool
       True if success
   
    Notes
    -----
    Success is assumed as a SHELX CC score of >= SHELXSUCCESS

    """
    directory, script = os.path.split(script_path)
    scriptname = os.path.splitext(script)[0]
    rfile = os.path.join(directory, 'search_' + scriptname + '_mrbump', 'results', 'resultsTable.pkl')
    if os.path.isfile(rfile):
        results = ResultsSummary().processMrbumpPkl(rfile)
        best = ResultsSummary.sortResultsStatic(results)[0]
        return ResultsSummary.jobSucceeded(best)
    else:
        return False


def finalSummary(amoptd):
    """Print a final summary of the job"""
    
    mrbump_data = amoptd['mrbump_results']
    if not mrbump_data:
        return "Could not find any MRBUMP results in directory: {0}!".format(amoptd['mrbump_dir'])
    
    if 'ensembles_data' in amoptd and not (amoptd['ideal_helices'] or amoptd['homologs'] or amoptd['single_model_mode']):
        results = []
        # Merge dictionaries together
        ensembles_data = amoptd['ensembles_data']
        for mrb in mrbump_data:
            d = copy.copy(mrb)
            for ed in ensembles_data:
                if ed['name'] == d['ensemble_name']:
                    d.update(ed)
                    results.append(d)
        keys = ['ensemble_name', 'Solution_Type', 'MR_program']
        keys += _resultsKeys(results)
        keys += ['subcluster_num_models', 'num_residues']
    else:
        results = mrbump_data
        keys = ['name', 'Solution_Type', 'MR_program']
        keys += _resultsKeys(results)
        
    resultsTable = []
    resultsTable.append(keys)
    for result in results:
        resultLine = []
        for k in keys:
            resultLine.append(result[k])
        resultsTable.append(resultLine)

    # Format the results
    table = printTable.Table()
    summary = table.pprint_table(resultsTable)

    r = "\n\nOverall Summary:\n\n"
    r += summary
    if len(results) and "MR_directory" in results[0]:
        r += '\nBest Molecular Replacement results so far are in:\n\n'
        r += str(results[0]["MR_directory"])
        r += '\n\n'
    return r


def job_unfinished(job_dict):
    if not 'Solution_Type' in job_dict: return True
    return job_dict['Solution_Type'] == "unfinished" or job_dict['Solution_Type'] == "no_job_directory"


def purge_MRBUMP(amoptd):
    """Remove as much as possible from a MRBUMP directory whilst keeping valid results"""
    mkey = 'mrbump_dir'
    if mkey not in amoptd or amoptd[mkey] is None:
        return
    mrbump_dir = amoptd[mkey]
    suffixes = ['.pdb', '.mtz', '.log', '.sh', '.mrbump']
    if os.path.isdir(mrbump_dir):
        for f in os.listdir(mrbump_dir):
            _, suffix = os.path.splitext(f)
            if suffix in suffixes:
                os.remove(os.path.join(mrbump_dir, f))
    return


def set_success_criteria(amoptd):
    """Set the module-level success criteria from an AMPLE job dictionary"""
    for criteria in ['SHELXE_CC', 'SHELXE_ACL']:
        amopt_prefix = 'early_terminate_'
        module_prefix = 'SUCCESS_'
        amopt_key = amopt_prefix + criteria
        if amopt_key in amoptd and amoptd[amopt_key] is not None:
            module_criteria = module_prefix + criteria
            logger.debug('Updating MRBUMP success criteria \'%s\' to: %s', module_criteria, amoptd[amopt_key])
            globals()[module_criteria] = amoptd[amopt_key]
        

def unfinished_scripts(amoptd):
    """See if there are any unfinished mrbump jobs in a mrbump directory and return a list of the scripts"""
    
    if not 'mrbump_dir' in amoptd or amoptd['mrbump_dir'] is None or not os.path.isdir(amoptd['mrbump_dir']): return []

    amoptd['mrbump_results'] = ResultsSummary().extractResults(amoptd['mrbump_dir'])
    if not len(amoptd['mrbump_results']): return []
    
    scripts = []
    for r in [ r for r in amoptd['mrbump_results'] if job_unfinished(r) ]:
        scripts.append(os.path.join(amoptd['mrbump_dir'], r['ensemble_name']+ample_util.SCRIPT_EXT))
    return scripts


def write_mrbump_files(ensemble_pdbs, amoptd, job_time=MRBUMP_RUNTIME, ensemble_options=None, directory=None):
    """Write the MRBUMP job files for all the ensembles.

    Arguments:
    ensemble_pdbs -- list of the ensembles, each a single pdb file.
    amoptd -- dictionary with job options.
    job_time -- maximum permissible runtime (mainly used for batch queueing systems).
    ensemble_options -- dictionary with ensemble-specific keywords e.g. ensemble_options[ensemble_name] = {'ncopies' : ncopies}
    directory -- working directory to write files to.
    """
    if not directory: directory = os.getcwd()
    
    job_scripts = []
    keyword_options = {}
    for ensemble_pdb in ensemble_pdbs:
        name = os.path.splitext(os.path.basename(ensemble_pdb))[0] # Get name from pdb path
        
        # Get any options specific to this ensemble
        if ensemble_options and name in ensemble_options: keyword_options = ensemble_options[name]
        
        # Generate dictionary with all the options for this job and write to keyword file
        keyword_dict = mrbump_cmd.keyword_dict(ensemble_pdb, name, amoptd, keyword_options)
        keyword_file = os.path.join(directory,name+'.mrbump')
        keyword_str = mrbump_cmd.mrbump_keyword_file(keyword_dict)
        with open(keyword_file,'w') as f: f.write(keyword_str)
        
        script = write_jobscript(name,
                                 keyword_file,
                                 amoptd,
                                 directory = directory,
                                 job_time = job_time)
        job_scripts.append(script)
            
    if not len(job_scripts):
        raise RuntimeError("No job scripts created!")
    
    return job_scripts


def write_jobscript(name, keyword_file, amoptd, directory=None, job_time=86400, extra_options={}):
    """
    Create the script to run MrBump for this PDB.
    """
    if not directory: directory = os.getcwd()
        
    # Next the script to run mrbump
    script_path = os.path.abspath(os.path.join(directory, name+ample_util.SCRIPT_EXT))
    with open(script_path, "w") as job_script:
        # Header
        if not sys.platform.startswith("win"):
            script_header = '#!/bin/sh\n'
            script_header += '[[ ! -d $CCP4_SCR ]] && mkdir $CCP4_SCR\n\n'
            job_script.write(script_header)
        
        # Get the mrbump command-line
        jobcmd = mrbump_cmd.mrbump_cmd(name, amoptd['mtz'], amoptd['mr_sequence'], keyword_file)
        job_script.write(jobcmd)
        
    # Make executable
    os.chmod(script_path, 0o777)
    logger.debug("Wrote MRBUMP script: {0}".format(script_path))
    return script_path

if __name__ == "__main__":
    if not len(sys.argv) == 2: 
        mrbump_dir = os.getcwd()
    else:
        mrbump_dir = os.path.join(os.getcwd(), sys.argv[1])
        
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    r = ResultsSummary()
    print(r.summariseResults(mrbump_dir))
