#!/usr/bin/env ccp4-python

import copy
import cPickle
import glob
import logging
import os
import shutil
import sys

# Hack to make sure we can find the modules we need
if __name__ == "__main__":
    root = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2])
    sys.path.insert(0, os.path.join(root, "scripts"))

# Our imports
from ample.util import ample_util
from ample.util import mrbump_cmd
from ample.util import printTable

# MRBUMP imports
if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
mrbumpd = os.path.join(os.environ['CCP4'], "share", "mrbump", "include", "parsers")
#mrbumpd = "/opt/mrbump-trunk/include/parsers"
sys.path.insert(0, mrbumpd)
import parse_arpwarp
import parse_buccaneer
import parse_phaser

TOP_KEEP = 3 # How many of the top shelxe/phaser results to keep for the gui
MRBUMP_RUNTIME = 172800 # allow 48 hours for each mrbump job

# We need a null logger so that we can be used without requiring a logger
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

LOGGER = logging.getLogger(__name__)

class ResultsSummary(object):
    """
    Summarise the results for a series of MRBUMP runs
    """

    def __init__(self):
        self.results = []
        # Add Null logger so we can be used without requiring a logger
        LOGGER.addHandler(NullHandler())
        self.pname = "archive"
        self.pdir = None
        self.success = False
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
            
        #
        # SHELXE PROCESSING
        #
        #
        # Buccaneer Rebuild Processing
        #
        buccaneerLog = os.path.join(mrDir,
                                     "build/shelxe/rebuild/buccaneer",
                                     "buccaneer.log")
    
        bp = parse_buccaneer.BuccaneerLogParser()
        if os.path.isfile(buccaneerLog):
            bp.parse(buccaneerLog)
            result["SXRBUCC_final_Rfree"] = bp.finalRfree
            result["SXRBUCC_final_Rfact"] = bp.finalRfact

        #
        # Arpwarp Rebuild Processing
        #
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
        # END
        
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
    
    def _extractOld(self, mrbump_dir):
        """Recreate a list of the jobs that have been purged"""
        old_results = {}
        self.pdir = os.path.join(mrbump_dir, self.pname)
        if not os.path.isdir(self.pdir): os.mkdir(self.pdir)
        pkls = glob.glob(os.path.join(self.pdir, "*.pkl"))
        if pkls:
            for p in pkls:
                with open(p) as f: d = cPickle.load(f)
                old_results[d['ensemble_name']] = d   
        return old_results

    def extractResults(self, mrbump_dir, purge=False):
        if not mrbump_dir or not os.path.isdir(mrbump_dir):
            raise RuntimeError,"Cannot find mrbump_dir: {0}".format(mrbump_dir)
        old_results = {}
        if purge: old_results = self._extractOld(mrbump_dir)
        results = self._extractResults(mrbump_dir, archived_ensembles=old_results.keys())
        
        if purge:
            self._purgeFailed(results)
            results += old_results.values()
            
        results = self.sortResults(results)
        
        self.success = any([jobSucceeded(r) for r in results])
        self.results = results
        return self.results

    def _extractResults(self, mrbump_dir, archived_ensembles=None):
        """
        Find the results from running MRBUMP and sort them
        """
        mrbump_dir = os.path.abspath(mrbump_dir)
        if not os.path.isdir(mrbump_dir):
            LOGGER.warn("extractResults - is not a valid directory: {0}".format(mrbump_dir))
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
            LOGGER.warn("Could not extract any results from directory: {0}".format(mrbump_dir))
            return []
        
        # reset any results
        results = []
        failed = {}  # dict mapping failures to what went wrong - need to process at the end
        for ensemble in ensembles:
            # Skip ones that we've archived
            if archived_ensembles and ensemble in archived_ensembles: continue
            
            # Check job directory
            jobDir = os.path.join(mrbump_dir, 'search_' + ensemble + '_mrbump')
            if not os.path.isdir(jobDir): jobDir = os.path.join( mrbump_dir, 'search_'+ensemble )
            if not os.path.isdir(jobDir):
                # As we call this every time we monitor a job running, we don't want to print this out all the time
                # LOGGER.debug("Missing job directory: {0}".format(jobDir))
                failed[ ensemble ] = "no_job_directory"
                continue

            LOGGER.debug(" -- checking directory for results: {0}".format(jobDir))
            # Check if finished
            if not os.path.exists(os.path.join(jobDir, "results", "finished.txt")):
                LOGGER.debug("Found unfinished job: {0}".format(jobDir))
                failed[ ensemble ] = "unfinished"
                continue

            # Check resultsTable.dat
            resultsDict = os.path.join(jobDir, "results", "resultsTable.pkl")
            resultsTable = os.path.join(jobDir, "results", "resultsTable.dat")
            if os.path.isfile(resultsDict):
                results += self.processResultsPkl(resultsDict)
            elif os.path.isfile(resultsTable):
                results += self.parseTableDat(resultsTable)
            else:
                LOGGER.debug(" -- Could not find results files: {0} or {1}".format(resultsDict, resultsTable))
                failed[ ensemble ] = "missing-results-file"
                continue

        # Process the failed results
        if failed: results += self._processFailed(mrbump_dir, failed)
        if not len(results): LOGGER.warn("Could not extract any results from directory: {0}".format(mrbump_dir))
        return results
    
    def parseTableDat(self, tfile):
        """Read a resultsTable file and return a list of MrBump results objects"""

        # Extract the various components from the path
        tfile = os.path.abspath(tfile)
        paths = tfile.split(os.sep)
        jobDir = os.path.abspath(os.sep.join(paths[:-2]))
        if paths[-3].endswith('_mrbump'): ensemble = paths[-3][7:-7]  #  remove search_..._mrbump: 'search_All_atom_trunc_0.005734_rad_1_mrbump'
        else: ensemble = paths[-3][7:]

        # List of all the possible column titles and their result object attributes
        title2key = {
                     'Model_Name'       : 'name',
                     'MR_Program'       : 'MR_program',
                     'Solution_Type'    : 'Solution_Type',
                     'final_Rfact'      : 'final_Rfact',
                     'final_Rfree'      : 'final_Rfree',
                     'SHELXE_CC'        : 'SHELXE_CC',
                     'SHELXE_ACL'       : 'SHELXE_ACL',
                     'Bucc_final_Rfact' : 'BUCC_final_Rfact',
                     'Bucc_final_Rfree' : 'BUCC_final_Rfree',
                     'ARP_final_Rfact'  : 'ARP_final_Rfact',
                     'ARP_final_Rfree'  : 'ARP_final_Rfree',
                      }

        results = []
        header = None
        nfields = None
        # Read results table to get the results
        for line in open(tfile):

            # Create a result object for each line in the output file
            result = self.createDict()
            result['ensemble_name'] = ensemble
            result['Search_directory'] = jobDir

            line = line.strip()
            if not header:
                # Processing header
                header = line.split()
                nfields = len(header)  # count as check
                for f in header:
                    # Map the data fields to their titles
                    if f not in title2key.keys():
                        LOGGER.critical("jobDir {0}: Problem with field {1} in headerline: {2}".format(jobDir, f, line))
                        result['Solution_Type'] = "problem-header-file.dat"
                        self._getUnfinishedResult(result)
                        results.append(result)
                        return results
                continue
                # End header processing

            fields = line.split()
            if len(fields) != nfields:
                msg = "jobDir {0}: Problem getting dataline: {1}".format(jobDir, line)
                LOGGER.debug(msg)
                result['Solution_Type'] = "corrupted-data-tfile.dat"
                self._getUnfinishedResult(result)
                results.append(result)
                continue

            # The headers tell us what attribute is in each column, so we use these and the header2attr dict to
            # set the results
            for i, title in enumerate(header):
                v = fields[i]
                if v == '--': v = None  # non-valid values in table are indicated by --
                result[title2key[title]] = v

            dirName = result['name'][:-6]
            result["MR_directory"] = os.path.join(jobDir, 'data', dirName, 'unmod', 'mr', result["MR_program"].lower())
            
            # See which pdb files were created
            if result["MR_program"] == 'PHASER':
                phaserPdb = os.path.join(result["MR_directory"],
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result["MR_program"].lower(), result['ensemble_name']))
                if os.path.isfile(phaserPdb):
                    result['PHASER_pdbout'] = phaserPdb
            elif result.program == 'molrep':
                molrepPdb = os.path.join(result["MR_directory"],
                                       "refine",
                                       "{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result["MR_program"].lower(), result['ensemble_name']))
                if os.path.isfile(molrepPdb):
                    result['MOLREP_pdbout'] = molrepPdb

            refmacPdb = os.path.join(result["MR_directory"],
                                   'refine',
                                   "refmac_" + result["MR_program"].lower() + "_loc0_ALL_" + result['ensemble_name'] + "_UNMOD.pdb")
            if os.path.isfile(refmacPdb):
                result["REFMAC_pdbout"] = refmacPdb
                
            shelxePdb = os.path.join(result["MR_directory"], 'build', 'shelxe',
                                   "shelxe_" + result["MR_program"].lower() + "_loc0_ALL_" + result['ensemble_name'] + "_UNMOD.pdb")
            if os.path.isfile(shelxePdb):
                result["SHELXE_pdbout"] = shelxePdb
                
            buccaneerPdb = os.path.join(result["MR_directory"], 'build', 'shelxe', 'rebuild', 'buccaneer',
                                      "buccSX_output.pdb")
            if os.path.isfile(buccaneerPdb):
                result["SXRBUCC_pdbout"] = buccaneerPdb
                
            arpWarpPdb = os.path.join(result["MR_directory"], 'build', 'shelxe', 'rebuild', 'arpwarp',
                                      "refmacSX_output_warpNtrace.pdb")
            if os.path.isfile(arpWarpPdb):
                result["SXRARP_pdbout"] = arpWarpPdb
            
            self.analyseResult(result)
            results.append(result)

        return results

    def processResultsPkl(self, resultsPkl):
        """Process dictionary
        """
        with open(resultsPkl) as f:
            rD = cPickle.load(f)
        if not rD: return []
        
        results = []
        for name, d1 in rD.iteritems():
            for mrprog, d2 in d1.iteritems():
                # Check if all the entries are None - means job didn't run.
                # Should probably think of a better way to spot that (Search_directory is always set)
                if not any([v for k, v in d2.iteritems() if k != 'Search_directory']): continue
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
            
        LOGGER.debug("Added {0} MRBUMP result failures".format(len(failed)))
        return results
    
    def _purgeFailed(self, results):
        """Remove any jobs that don't pass the keep criteria and archive their job dictionaries"""
        # Skip any that are unfinished
        completed = [ r for r in results if not(job_unfinished(r)) ]
        if completed:
            # Keep the top TOP_KEEP SHELXE_CC and PHASER_TFZ - these could be the same jobs and we may not even
            # have TOP_KEEP completed
            to_keep = []
            for r in self.sortResults(completed, prioritise='SHELXE_CC')[0:min(len(completed),TOP_KEEP)]:
                if r not in to_keep: to_keep.append(r)
            for r in self.sortResults(completed, prioritise='PHASER_TFZ')[0:min(len(completed),TOP_KEEP)]:
                if r not in to_keep: to_keep.append(r)

            # Remove the directories and archive the dictionaries
            for r in completed:
                if r not in to_keep:
                    pkl = os.path.join(self.pdir, "{0}.pkl".format(r['ensemble_name']))
                    with open(pkl, 'w') as f: cPickle.dump(r, f)
                    shutil.rmtree(r['Search_directory'])
        
#         for r in results:
#             if r['Solution_Type'] == "unfinished" or r['Solution_Type'] == "no_job_directory": continue
#             to_keep.append(r)
#             if not (jobSucceeded(r) or r['Solution_Type'] is "MARGINAL") :
#                 pkl = os.path.join(self.pdir, "{0}.pkl".format(r['ensemble_name']))
#                 with open(pkl, 'w') as f: cPickle.dump(r, f)
#                 shutil.rmtree(r['Search_directory'])
        return

    def results_table(self, results):
        resultsTable = []
        keys = ['ensemble_name', 'MR_program', 'Solution_Type']
        keys += _resultsKeys(results)
        
        resultsTable.append(keys)
        for r in results: resultsTable.append([r[k] for k in keys])
        return resultsTable

    def sortResults(self, results, prioritise=None):
        """
        Sort the results
        """
        # Check each result to see what attributes are set and use this to work out how we rate this run
        
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
        if SHELXE and not prioritise == "PHASER_TFZ":
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
            # Now sort by the key
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
            if p in top.keys(): k = p
        assert k,"Missing search directory key in results dictionary"
        if top[k]:
            r += '\nBest Molecular Replacement results so far are in:\n\n'
            r += self.results[0]["MR_directory"]
        r += '\n\n'

        return r

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
    rfile = os.path.join(directory, 'search_' + scriptname + '_mrbump', 
                         'results', 'resultsTable.pkl')
    # print "{0} checking for file: {1}".format(multiprocessing.current_process().name,rfile)
    if not os.path.isfile(rfile):
        # print "{0} cannot find results file: {1}".format(multiprocessing.current_process().name,rfile)
        return False
    
    # Results summary object to parse table file
    mrbR = ResultsSummary()
    
    # Put into order and take top one
    results = mrbR.processResultsPkl(rfile)
    mrbR.sortResults(results)
    return jobSucceeded(results[0])

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

def jobSucceeded(job_dict):
    PHASER_TFZ = 8.0
    PHASER_LLG = 120
    RFREE = 0.4
    SHELXE_CC = 25.0
    SHELXE_ACL = 10
    success = False
    if 'SHELXE_CC' in job_dict and job_dict['SHELXE_CC'] and float(job_dict['SHELXE_CC']) >= SHELXE_CC and \
       'SHELXE_ACL' in job_dict and job_dict['SHELXE_ACL'] and float(job_dict['SHELXE_ACL']) >= SHELXE_ACL:
        success = True
    elif 'BUCC_final_Rfact' in job_dict and job_dict['BUCC_final_Rfact'] and float(job_dict['BUCC_final_Rfact']) <= RFREE:
        success = True
    elif 'ARP_final_Rfree' in job_dict and job_dict['ARP_final_Rfree'] and float(job_dict['ARP_final_Rfree']) <= RFREE:
        success = True
    elif 'REFMAC_Rfree' in job_dict and job_dict['REFMAC_Rfree'] and float(job_dict['REFMAC_Rfree']) <= RFREE:
        success = True
    elif 'PHASER_LLG' in job_dict and 'PHASER_TFZ' in job_dict and job_dict['PHASER_LLG'] and job_dict['PHASER_TFZ'] and \
    float(job_dict['PHASER_LLG']) >= PHASER_LLG and float(job_dict['PHASER_TFZ']) >= PHASER_TFZ:
        success = True
    return success

def job_unfinished(job_dict):
    if not 'Solution_Type' in job_dict: return True
    return job_dict['Solution_Type'] == "unfinished" or job_dict['Solution_Type'] == "no_job_directory"

def unfinished_scripts(amoptd):
    """See if there are any unfinished mrbump jobs in a mrbump directory and return a list of the scripts"""
    
    if not 'mrbump_dir' in amoptd or amoptd['mrbump_dir'] is None or not os.path.isdir(amoptd['mrbump_dir']): return []

    amoptd['mrbump_results'] = ResultsSummary().extractResults(amoptd['mrbump_dir'])
    if not len(amoptd['mrbump_results']): return []
    
    scripts = []
    for r in [ r for r in amoptd['mrbump_results'] if job_unfinished(r) ]:
        #print "DIR ", r['MR_directory']
        #print "DIR2 ", r['Search_directory']
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
        msg = "No job scripts created!"
        logging.critical(msg)
        raise RuntimeError, msg
    
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
    LOGGER.debug("Wrote MRBUMP script: {0}".format(script_path))

    return script_path


if __name__ == "__main__":
    if not len(sys.argv) == 2: 
        mrbump_dir = os.getcwd()
    else:
        mrbump_dir = os.path.join(os.getcwd(), sys.argv[1])
        
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    r = ResultsSummary()
    print r.summariseResults(mrbump_dir)

