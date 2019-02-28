#! /usr/bin/env ccp4-python

import argparse
import os
import pickle
import sys

from mrbump.cluster import refmacEXE
from mrbump.building import MRBUMP_Shelxe
from mrbump.building import MRBUMP_Buccaneer
from mrbump.building import MRBUMP_ARPwARP
from mrbump.parsers import parse_buccaneer
from mrbump.parsers import parse_arpwarp
from mrbump.tools import MRBUMP_phs2mtz, make_dictionary


from mrbump.structures.Matches import Match_struct, SEQ_match, Process_Score
from mrbump.structures.Model_struct import MR_setup
from mrbump.file_info.MRBUMP_target_info import TargetInfo
from mrbump.modelling.Model import Modelling
from mrbump.tools.make_dictionary import makeSearchDict
from mrbump.initialisation.MRBUMP_initialise import Initialise



class ShelxePlusRebuild:
    """Run shelxe and rebuild"""

    def __init__(self):
        self.logfile = ""
        self.DB_fail = False
        if os.name == "nt":
            self.pdbcurEXE = os.path.join(os.environ["CCP4_BIN"], "pdbcur.exe")
            self.refmacEXE = os.path.join(os.environ["CCP4_BIN"], "refmac5.exe")
        else:
            self.pdbcurEXE = os.path.join(os.environ["CCP4_BIN"], "pdbcur")
            self.refmacEXE = os.path.join(os.environ["CCP4_BIN"], "refmac5")
        try:
            self.debug = eval(os.environ['MRBUMP_DEBUG'])
        except:
            self.debug = False

    def shelxe_trace(self, init, mstat, model, target_info):

        shelxe_job = MRBUMP_Shelxe.Shelxe() 
        #shelxe_dir = os.path.join(model.model_directory, "mr", "phaser", "build", "shelxe")
        shelxe_dir = os.path.join(model.model_directory, "shelxe")
        shelxe_job.setShelxWorkingDIR(shelxe_dir)
        shelxe_job.setShelxeLogFile(os.path.join(shelxe_dir, "shelxe_run.log"))
        shelxe_PDBINfile = model.refmac_phaser_PDBfile
        shelxe_MTZINfile = model.refmac_phaser_MTZINfile
        # Run SHELXE to do phase improvement and a c-alpha trace
        # Build up list of argument and keyword arguments
        shelxe_args = [target_info.solvent, target_info.resolution, "PHASER", shelxe_PDBINfile, shelxe_MTZINfile]
        shelxe_kwargs = {
            'pdboutFile': model.shelxe_phaser_PDBfile,
            'phsoutFile': model.shelxe_phaser_PHSfile,
            'traceCycles': init.keywords.SCYCLES,
            'fp': init.workingData.mtz_coldata['F'],
            'sigfp': init.workingData.mtz_coldata['SIGF'],
            'free': init.workingData.mtz_coldata['FREE'],
            'shelxeEXE': init.keywords.SHLXEXE
        }
        if not init.keywords.PDBNATIVE is None:
            shelxe_kwargs['native'] = init.keywords.PDBNATIVE

        shelxe_job.runShelxe(*shelxe_args, **shelxe_kwargs)

        # Convert the phs file to an mtz
        if os.path.isfile(model.shelxe_phaser_PHSfile) and os.path.isfile(model.shelxe_phaser_PDBfile):
            p2m = MRBUMP_phs2mtz.PHS2MTZ()
            p2m.phs2mtz(
                model.shelxe_phaser_PHSfile,
                model.shelxe_phaser_PDBfile,
                model.shelxe_phaser_MTZfile,
                shelxe_job.workingDIR,
                hklref=init.workingData.hklin,
                f=init.workingData.mtz_coldata['F'],
                sigf=init.workingData.mtz_coldata['SIGF'],
                freeLabel=init.workingData.mtz_coldata['FREE'],
                resolution=target_info.resolution,
                debug=self.debug)
        else:
            sys.stdout.write("%s: warning: no PHS file found for SHELXE (after phaser):\n  %s" %
                             (model.name, model.shelxe_phaser_PHSfile))

        # Set the SHELXE CC value and average chain length for partial model
        shelxe_job.parseLog()
        cc = shelxe_job.resultsDict['SHELXE_CC']
        acl = shelxe_job.resultsDict['SHELXE_ACL']
        if cc and acl:  # Only set if found as are initialised to 0.0
            model.shelxe_phaser_CCscore = cc
            model.shelxe_phaser_AvgChainLen = acl
        mstat.results_dict[model.name]["PHASER"]["SHELXE_CC"] = cc
        mstat.results_dict[model.name]["PHASER"]["SHELXE_ACL"] = acl
        mstat.results_dict[model.name]["PHASER"]["SHELXE_MCL"] = shelxe_job.resultsDict['SHELXE_MCL']
        mstat.results_dict[model.name]["PHASER"]["SHELXE_NC"] = shelxe_job.resultsDict['SHELXE_NC']
        mstat.results_dict[model.name]["PHASER"]["SHELXE_time"] = shelxe_job.resultsDict['SHELXE_time']
        mstat.results_dict[model.name]["PHASER"]["SHELXE_version"] = shelxe_job.resultsDict['SHELXE_version']
        mstat.results_dict[model.name]["PHASER"]["SHELXE_wMPE"] = shelxe_job.resultsDict['SHELXE_wMPE']
        mstat.results_dict[model.name]["PHASER"]["SHELXE_os"] = shelxe_job.resultsDict['SHELXE_os']

        #sys.stdout.write("%s: warning: no log file found for SHELXE (after phaser):\n  %s" % (model.name, model.shelxe_phaser_logfile))

        return shelxe_job

    def rebuild_buccaneer(self, init, mstat, model, target_info, rebuild_dir, refmacSX_MTZout, refmacSX_PDBout):
        # Run Buccaneer to build on the refined SHELXE trace
        buccaneerSX = MRBUMP_Buccaneer.Buccaneer()
        os.mkdir(os.path.join(rebuild_dir, "buccaneer"))  # Set the output MTZ and PDB file from Buccaneer
        model.buccSHELXE_phaser_PDBfile = os.path.join(rebuild_dir, "buccaneer", "buccSX_output.pdb")
        model.buccSHELXE_phaser_MTZfile = os.path.join(rebuild_dir, "buccaneer", "buccaneer_pipeline", "refine.mtz")
        buccSHELXE_phaser_logfile = os.path.join(rebuild_dir, "buccaneer", "buccaneer.log")
        # Run Buccaneer
        if init.keywords.USEPHS:
            buccaneerSX.runBuccaneer(
                init.seqin,
                model.shelxe_phaser_MTZfile,
                model.buccSHELXE_phaser_PDBfile,
                os.path.join(rebuild_dir, "buccaneer"),
                'F',
                'SIGF',
                init.workingData.mtz_coldata['FREE'],
                "PHI_SHELXE",
                "FOM_SHELXE",
                cycles=5,
                pdbinFile=model.shelxe_phaser_PDBfile,
                WEB_PATH_START=init.WEB_PATH_START)
        else:
            buccaneerSX.runBuccaneer(
                init.seqin,
                refmacSX_MTZout,
                model.buccSHELXE_phaser_PDBfile,
                os.path.join(rebuild_dir, "buccaneer"),
                init.workingData.mtz_coldata['F'],
                init.workingData.mtz_coldata['SIGF'],
                init.workingData.mtz_coldata['FREE'],
                "PHIC",
                "FOM",
                "FWT",
                "PHWT",
                cycles=5,
                pdbinFile=refmacSX_PDBout,
                WEB_PATH_START=init.WEB_PATH_START)

        # Add the bits to the results dictionary
        mstat.results_dict[model.name]["PHASER"]["SXRBUCC_pdbout"] = model.buccSHELXE_phaser_PDBfile
        mstat.results_dict[model.name]["PHASER"]["SXRBUCC_mtzout"] = model.buccSHELXE_phaser_MTZfile
        if os.path.isfile(buccSHELXE_phaser_logfile):
            mstat.results_dict[model.name]["PHASER"]["SXRBUCC_logfile"] = buccSHELXE_phaser_logfile
            bp = parse_buccaneer.BuccaneerLogParser(buccSHELXE_phaser_logfile)
            mstat.results_dict[model.name]["PHASER"]["SXRBUCC_version"] = bp.version
            mstat.results_dict[model.name]["PHASER"]["SXRBUCC_final_Rfact"] = bp.finalRfact
            mstat.results_dict[model.name]["PHASER"]["SXRBUCC_final_Rfree"] = bp.finalRfree
        return

    def rebuild_arpwarp(self, init, mstat, model, target_info, rebuild_dir, refmacSX_MTZout, refmacSX_PDBout):
        # Run ARP/wARP to build on the refined SHELXE trace
        arpwarpSX = MRBUMP_ARPwARP.Arpwarp()
        os.mkdir(os.path.join(rebuild_dir, "arpwarp"))
        # Set the output MTZ and PDB file from ARP/wARP
        if init.keywords.USEPHS:
            model.arpSHELXE_phaser_PDBfile = os.path.join(
                rebuild_dir, "arpwarp",
                "%s_warpNtrace.pdb" % os.path.splitext(os.path.split(model.shelxe_phaser_MTZfile)[1])[0])
            model.arpSHELXE_phaser_MTZfile = os.path.join(
                rebuild_dir, "arpwarp",
                "%s_warpNtrace.mtz" % os.path.splitext(os.path.split(model.shelxe_phaser_MTZfile)[1])[0])
        else:
            model.arpSHELXE_phaser_PDBfile = os.path.join(rebuild_dir, "arpwarp", "refmacSX_output_warpNtrace.pdb")
            model.arpSHELXE_phaser_MTZfile = os.path.join(rebuild_dir, "arpwarp", "refmacSX_output_warpNtrace.mtz")
        arpSHELXE_phaser_logfile = os.path.join(rebuild_dir, "arpwarp", "arpwarp.log")
        # Set the number of residues in the asu
        noResiduesASU = target_info.no_of_res * target_info.no_of_mols
        # Run ARP/wARP
        if init.keywords.USEPHS:
            arpwarpSX.runARPwARP(
                init.seqin,
                model.shelxe_phaser_MTZfile,
                os.path.join(rebuild_dir, "arpwarp"),
                'F',
                'SIGF',
                init.workingData.mtz_coldata['FREE'],
                "PHI_SHELXE",
                "FOM_SHELXE",
                noResiduesASU,
                cycles=5,
                pdbinFile=model.shelxe_phaser_PDBfile)
        else:
            arpwarpSX.runARPwARP(
                init.seqin,
                refmacSX_MTZout,
                os.path.join(rebuild_dir, "arpwarp"),
                init.workingData.mtz_coldata['F'],
                init.workingData.mtz_coldata['SIGF'],
                init.workingData.mtz_coldata['FREE'],
                "PHIC",
                "FOM",
                noResiduesASU,
                cycles=5,
                pdbinFile=refmacSX_PDBout)
        # Add the bits to the results dictionary
        mstat.results_dict[model.name]["PHASER"]["SXRARP_pdbout"] = model.arpSHELXE_phaser_PDBfile
        mstat.results_dict[model.name]["PHASER"]["SXRARP_mtzout"] = model.arpSHELXE_phaser_MTZfile
        if os.path.isfile(arpSHELXE_phaser_logfile):
            mstat.results_dict[model.name]["PHASER"]["SXRARP_logfile"] = arpSHELXE_phaser_logfile
            ap = parse_arpwarp.ArpwarpLogParser(arpSHELXE_phaser_logfile)
            mstat.results_dict[model.name]["PHASER"]["SXRARP_version"] = ap.version
            mstat.results_dict[model.name]["PHASER"]["SXRARP_final_Rfact"] = ap.finalRfact
            mstat.results_dict[model.name]["PHASER"]["SXRARP_final_Rfree"] = ap.finalRfree
        return

    def submit_serial(self, init, mstat, model, target_info):
        #=================================================================================#
        #================== Run MR and refinement on the local machine ===================#
        #=================================================================================#

        ###############################################################################
        ############################# Running Phaser ##################################
        ###############################################################################
        self.logfile = os.path.join(model.model_directory, "phaser_refmac_job.log")

        # C-alpha tracing with Shelxe
        shelxe_job = self.shelxe_trace(init, mstat, model, target_info)
        if not init.keywords.SXREBUILD:
            return

        # Rebuilding the shelxe trace
        if not os.path.isfile(model.shelxe_phaser_PDBfile) or not os.path.isfile(model.shelxe_phaser_PHSfile):
            sys.stdout.write("Warning: No output PDB from SHELXE, can't run SHELXE rebuild step\n")
            sys.stdout.write("\n")
            return

        if not shelxe_job.succeeded():
            sys.stdout.write("Warning: SHELXE trace didn't pass success criteria so not continuing with rebuilding\n")
            sys.stdout.write("\n")
            return

        # Directory where all rebuilding is carried out
        rebuild_dir = os.path.join(shelxe_job.workingDIR, "rebuild")
        os.mkdir(rebuild_dir)

        # Take the SHELXE result and use it as a starting point in Buccaneer and/or ARP/wARP for building the target
        refmacSX_MTZout = None
        refmacSX_PDBout = None
        if not init.keywords.USEPHS:
            # If we are using the model file from SHELXE for phase information we need to run refmac on it
            # Run a round of refinement and model building with Buccaneer post SHELXE
            os.mkdir(os.path.join(rebuild_dir, "refine"))
            # Setup a Refmac run
            refmacSX_job = refmacEXE.RefmacEXE()
            # Setup the Refmac keywords
            refmacSX_job.add_keyword("LABIN FP=" + init.workingData.mtz_coldata['F'] + " SIGFP=" + init.workingData.
                                     mtz_coldata['SIGF'] + " FREE=" + init.workingData.mtz_coldata['FREE'])
            refmacSX_job.add_keyword('ridg DIST SIGM 0.02')
            refmacSX_job.add_keyword('weight auto')
            refmacSX_job.add_keyword('ncyc 30')
            refmacSX_job.add_keyword('make hydr no')
            refmacSX_job.add_keyword('END')
            # Set the output and log files for Refmac
            refmacSX_MTZout = os.path.join(rebuild_dir, "refine", "refmacSX_output.mtz")
            refmacSX_PDBout = os.path.join(rebuild_dir, "rebuild", "refine", "refmacSX_output.pdb")
            refmacSX_Logfile = os.path.join(rebuild_dir, "rebuild", "refine", "refmacSX.log")
            # Run refmac Restrained Refinement with Jelly Body
            refmacSX_job.run(
                model.refmac_phaser_MTZINfile,
                refmacSX_MTZout,
                model.shelxe_phaser_PDBfile,
                refmacSX_PDBout,
                refmacSX_Logfile,
                "Restrained Refinement with Jelly Body",
                refineDir=os.path.join(rebuild_dir, "refine"),
                script=os.path.join(rebuild_dir, "refine", "refmacSX_run_script.sh"),
                LITE=False,
                WEB_PATH_START=init.WEB_PATH_START)

        if init.keywords.SXRBUCC:
            self.rebuild_buccaneer(init, mstat, model, target_info, rebuild_dir, refmacSX_MTZout, refmacSX_PDBout)

        if init.keywords.SXRARPW:
            self.rebuild_arpwarp(init, mstat, model, target_info, rebuild_dir, refmacSX_MTZout, refmacSX_PDBout)
        return
 
def setup(search_dir):
    mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
    mrbump_incl = os.path.join(mrbump, "include")
    mrbump_data = os.path.join(mrbump, "data")
     
    os.environ['MRBUMP_DEBUG'] = "True"
    os.environ['MRBUMP_CLUSTER'] = "False"
    os.environ['PYTHONUNBUFFERED'] = "True"
    os.environ["GFORTRAN_UNBUFFERED_ALL"] = "Y"
    os.environ["MRBUMP_CLUSTBIN"] = os.path.join(mrbump_incl, "cluster")
    if 'CCP4_BIN' not in os.environ:
        os.environ['CCP4_BIN'] = os.path.join(os.environ['CCP4'], "bin")
    if os.path.isdir(os.path.join(os.environ['CCP4'], 'libexec')) and os.name=="nt":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + os.path.join(os.environ["CCP4"], "libexec")

    init = Initialise()
    init.setMRBUMP(mrbump)
    init.setMRBUMP_INCL(mrbump_incl)
    init.setMRBUMP_DATA(mrbump_data) 
    init.search_dir = os.path.abspath(search_dir)
    init.search_dict = makeSearchDict(search_dir)
    init.make_job_dir()
 
    mstat = Match_struct()
    mstat.search_dir = init.search_dir
    mstat.results_dir = init.results_dir
 
    # Need to make sure we don't try and write out an I2 report
    def mock_makeI2report(*args):
        pass
    mstat.makeI2report = mock_makeI2report

    return init, mstat


parser = argparse.ArgumentParser()
parser.add_argument('--name')
parser.add_argument('--resolution', type=float)
parser.add_argument('--seqin')
parser.add_argument('--hklin')
parser.add_argument('--pdb_native')
parser.add_argument('--f_label')
parser.add_argument('--sigf_label')
parser.add_argument('--free_label')
parser.add_argument('--refmac_mtz')
parser.add_argument('--refmac_pdb')

args = parser.parse_args()
# Put everything from the command-line into local variables
locals().update(vars(args))

# Main
search_dir = os.path.join(os.getcwd(), name)
init, mstat = setup(search_dir)

# Global keywords
init.keywords.SXREBUILD = True
init.keywords.SXRBUCC = True

# Specify info on target
target_info = TargetInfo()
# target_info.seqfile = 'FIXME'
# target_info.sequence = 'FIXME'
# target_info.mol_weight = 7071.09
# target_info.no_of_mols = 1
# target_info.solvent =
target_info.resolution = resolution
 
init.workingData.mtz_coldata['F'] = f_label
init.workingData.mtz_coldata['SIGF'] = sigf_label
init.workingData.mtz_coldata['FREE'] = free_label
init.workingData.hklin = hklin
init.keywords.PDBNATIVE = pdb_native
init.seqin = seqin


model = Modelling()
model.name = name
model.model_directory = os.path.join(init.search_dir, model.name)
mstat.results_dict[model.name] = make_dictionary.makeDict(init.search_dir)
shelxe_dir = os.path.join(model.model_directory, "shelxe")
os.makedirs(shelxe_dir)

model.refmac_phaser_PDBfile = refmac_pdb
model.refmac_phaser_MTZINfile = refmac_mtz
model.shelxe_phaser_PDBfile = os.path.join(shelxe_dir  ,"shelxe_phaser_" + model.name + ".pdb")
model.shelxe_phaser_PHSfile = os.path.join(shelxe_dir  ,"shelxe_phaser_" + model.name + ".phs")
model.shelxe_phaser_MTZfile = os.path.join(shelxe_dir  ,"shelxe_phaser_" + model.name + ".mtz")

spr = ShelxePlusRebuild()
spr.submit_serial(init, mstat, model, target_info)

init.resultspkl = os.path.join(init.search_dir, "results", "resultsTable.pkl")
with open(os.path.join(init.resultspkl), "w") as f:
    pickle.dump(mstat.results_dict, f)

