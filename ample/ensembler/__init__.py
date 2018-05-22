"""The AMPLE ensembler sub-package

This sub-package contains the source to apply AMPLE's
cluster-and-truncate approach to protein structure models.
It contains a series of different modules to apply varying
cluster and truncation algorithms to those models.
"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "14 Nov 2016"
__version__ = "1.0"

import collections
import glob
import iotbx.pdb
import logging
import operator
import os
import shutil
import sys

from abinitio import AbinitioEnsembler
from constants import SPICKER_TM
from homologs import HomologEnsembler
from single_model import SingleModelEnsembler
from ample.util import ample_util
from ample.util import argparse_util
from ample.util import exit_util
from ample.util import pdb_edit
from ample.util import printTable

logger = logging.getLogger(__name__)


def cluster_script(amoptd, python_path="ccp4-python"):
    """Create the script for ensembling on a cluster

    Parameters
    ----------
    amoptd : dict
       An AMPLE option dictionary
    python_path : str, optional
       The path to the CCP4 python executable

    Returns
    -------
    str
       The path to the cluster script

    """
    # write out script
    work_dir = amoptd['work_dir']
    script_path = os.path.join(work_dir, "submit_ensemble.sh")
    with open(script_path, "w") as job_script:
        job_script.write(ample_util.SCRIPT_HEADER + os.linesep)
        job_script.write("export CCP4_SCR=${TMPDIR}" + os.linesep) #Added by Ronan after issues on CCP4online server
        job_script.write("ccp4-python -m ample.ensembler -restart_pkl {0}".format(amoptd['results_path']) + os.linesep)

    # Make executable
    os.chmod(script_path, 0o777)
    return script_path


def create_ensembles(amoptd):
    """Create the ensembles using the values in the amoptd dictionary

    Parameters
    ----------
    amoptd : dict
       An AMPLE option dictionary

    """
    # Create instance of the ensembler
    ensembler = ensembler_factory(amoptd)

    ############################################################################
    # For a single model we don't need to use glob
    if not (amoptd['single_model'] or amoptd['models']):
        msg = 'AMPLE ensembler needs either a single_model or a list of models'
        exit_util.exit_error(msg, sys.exc_info()[2])
        if amoptd['single_model'] and not os.path.isfile(amoptd['single_model']):
            msg = 'Cannot find single_model pdb: {0}'.format(amoptd['single_model'])
            exit_util.exit_error(msg, sys.exc_info()[2])
        elif amoptd['models'] and len(amoptd['models'] < 2):
            msg = 'Not enough models provided for ensembling - use single_model_mode instead'
            exit_util.exit_error(msg, sys.exc_info()[2])

    models = list([amoptd['single_model']]) if amoptd['single_model_mode'] else amoptd['models']

    # Run ensemble creation
    ensembles = ensembler.generate_ensembles_from_amoptd(models, amoptd)

    ############################################################################
    # Hack to pull out the data - need to update code to work with ensemble objects rather than dictionaries
    amoptd['ensembles'] = [e.pdb for e in ensembles]
    amoptd['ensembles_data'] = [e.__dict__ for e in ensembles]

    # We need to let the main process know that we have succeeded as this module could be run on a cluster node with no link
    # to the parent process, so we create a file here indicating that we got this far and didn't die from an exception
    with open(amoptd['ensemble_ok'], 'w') as f:
        f.write('ok\n')

    # Delete all intermediate files if we're purging
    if amoptd['purge']:
        shutil.rmtree(ensembler.work_dir)
    return


def ensembler_factory(amoptd):
    """Return an ensembler object for the required ensembles

    Parameters
    ----------
    amoptd : dict
       An AMPLE option dictionary

    Returns
    -------
    :obj:`HomologEnsembler`, :obj:`SingleModelEnsembler`, or :obj:`AbinitioEnsembler`
       Instantiated ensembler object

    """
    if amoptd['homologs']:
        ensembler_class = HomologEnsembler
    elif amoptd['single_model_mode']:
        ensembler_class = SingleModelEnsembler
    else:
        ensembler_class = AbinitioEnsembler

    # Paths to ensembling directories
    ensembles_directory = os.path.join(amoptd['work_dir'], 'ensembles')
    amoptd['ensembles_directory'] = ensembles_directory
    work_dir = os.path.join(amoptd['work_dir'], 'ensemble_workdir')

    return ensembler_class(
        ensembles_directory=ensembles_directory,
        ensemble_max_models=amoptd['ensemble_max_models'],
        nproc=amoptd['nproc'],
        work_dir=work_dir,
        # Executables
        fast_protein_cluster_exe=amoptd['fast_protein_cluster_exe'],
        gesamt_exe=amoptd['gesamt_exe'],
        lsqkab_exe=amoptd['lsqkab_exe'],
        maxcluster_exe=amoptd['maxcluster_exe'],
        scwrl_exe=amoptd['scwrl_exe'],
        spicker_exe=amoptd['spicker_exe'],
        theseus_exe=amoptd['theseus_exe'],
    )


def collate_cluster_data(ensembles_data):
    """Collate the cluster data into a dictionary structure

    Parameters
    ----------
    ensembles_data : list, tuple
       A list of ensembles' data dictionaries

    Returns
    -------
    A dictionary with the collated data

    """

    clusters = {}  # Loop through all ensemble data objects and build up a data tree
    cluster_method = None
    cluster_score_type = None
    truncation_method = None
    percent_truncation = None
    side_chain_treatments = []
    for e in ensembles_data:
        if not cluster_method:
            cluster_method = e['cluster_method']
            cluster_score_type = e['cluster_score_type']
            percent_truncation = e['truncation_percent']
            truncation_method = e['truncation_method']
            # num_clusters = e['num_clusters']
        cnum = e['cluster_num']
        if cnum not in clusters:
            clusters[cnum] = {}
            clusters[cnum]['cluster_centroid'] = e['cluster_centroid']
            clusters[cnum]['cluster_num_models'] = e['cluster_num_models']
            clusters[cnum]['tlevels'] = {}
        tlvl = e['truncation_level']
        if tlvl not in clusters[cnum]['tlevels']:
            clusters[cnum]['tlevels'][tlvl] = {}
            clusters[cnum]['tlevels'][tlvl]['truncation_variance'] = e['truncation_variance']
            clusters[cnum]['tlevels'][tlvl]['num_residues'] = e['num_residues']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'] = {}
        srt = e['subcluster_radius_threshold']
        if srt not in clusters[cnum]['tlevels'][tlvl]['radius_thresholds']:
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt] = {}
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['num_models'] = e['subcluster_num_models']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'] = {}
        sct = e['side_chain_treatment']
        if sct not in side_chain_treatments:
            side_chain_treatments.append(sct)
        if sct not in clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct']:
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct] = {}
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['name'] = e['name']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['num_atoms'] = e['ensemble_num_atoms']

    return {
        'clusters': clusters,
        'cluster_method': cluster_method,
        'cluster_score_type': cluster_score_type,
        'truncation_method': truncation_method,
        'percent_truncation': percent_truncation,
        'side_chain_treatments': side_chain_treatments,
    }


def cluster_table_data(clusters, cluster_num, side_chain_treatments, header=True):
    """Generate a table containing the cluster data

    Parameters
    ----------
    clusters : list, tuple
       A list of cluster dictionaries
    cluster_num : int
       The cluster number of interest
    side_chain_treatments : list, tuple
       A list of side chain treatments
    header: bool
      Return a header line

    Returns
    -------
    list
       A list of table data lines

    """
    # FIX TO IGNORE side_chain_treatments
    # tdata = [("Name", "Truncation Level", u"Variance Threshold (\u212B^2)", "No. Residues", u"Radius Threshold (\u212B)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
    if header:
        tdata = [("Name", "Cluster", "Truncation Level", "Variance Threshold (A^2)", "No. Residues",
                  "Radius Threshold (A)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
    else:
        tdata = []
    for tl in sorted(clusters[cluster_num]['tlevels']):
        tvar = clusters[cluster_num]['tlevels'][tl]['truncation_variance']
        nresidues = clusters[cluster_num]['tlevels'][tl]['num_residues']
        for i, rt in enumerate(sorted(clusters[cluster_num]['tlevels'][tl]['radius_thresholds'])):
            nmodels = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['num_models']
            side_chain_treatments = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'].keys()
            for j, sct in enumerate(side_chain_treatments):
                name = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['name']
                num_atoms = clusters[cluster_num]['tlevels'][tl]['radius_thresholds'][rt]['sct'][sct]['num_atoms']
                if i == 0 and j == 0:  # change of radius
                    tdata.append((name, cluster_num, tl, tvar, nresidues, rt, nmodels, num_atoms, sct))
                elif i > 0 and j == 0:  # change of side_chain
                    tdata.append((name, "", "", "", "", rt, nmodels, num_atoms, sct))
                else:
                    tdata.append((name, "", "", "", "", "", "", num_atoms, sct))
    return tdata


def ensemble_summary(ensembles_data):
    """Print a summary of the ensembling process

    Parameters
    ----------
    ensembles_data : list, tuple
       A list of ensembles' data dictionaries

    Returns
    -------
    str

    """

    tableFormat = printTable.Table()

    d = collate_cluster_data(ensembles_data)
    clusters = d['clusters']
    num_clusters = len(clusters)

    rstr = "\n"
    rstr += "Ensemble Results\n"
    rstr += "----------------\n\n"
    rstr += "Cluster method: {0}\n".format(d['cluster_method'])
    rstr += "Cluster score type: {0}\n".format(d['cluster_score_type'])
    rstr += "Number of clusters: {0}\n".format(num_clusters)
    rstr += "Truncation method: {0}\n".format(d['truncation_method'])
    rstr += "Percent truncation: {0}\n".format(d['percent_truncation'])
    rstr += "Side-chain treatments: {0}\n".format(d['side_chain_treatments'])

    for cluster_num in sorted(clusters.keys()):
        rstr += "\n"
        rstr += "Cluster {0}\n".format(cluster_num)
        rstr += "Number of models: {0}\n".format(clusters[cluster_num]['cluster_num_models'])
        rstr += "Cluster centroid: {0}\n".format(clusters[cluster_num]['cluster_centroid'])
        rstr += "\n"
        tdata = cluster_table_data(clusters, cluster_num, d['side_chain_treatments'])
        rstr += tableFormat.pprint_table(tdata)

    rstr += "\nGenerated {0} ensembles\n\n".format(len(ensembles_data))

    return rstr


def get_ensembler_timeout(optd, tm_timeout=3600*8):
    """Set how long the ensembling should run based on the type of job being run"""
    timeout = optd['ensembler_timeout']
    if optd['cluster_method'] == SPICKER_TM:
        if int(timeout) < tm_timeout:
            timeout = tm_timeout
    return timeout


def import_ensembles(amoptd):
    """Import ensembles using their file paths

    Parameters
    ----------
    amoptd : dict
       An AMPLE option dictionary

    Returns
    -------
    list
       A list of absolute files paths of the ensembles

    """
    if not pdb_edit.check_pdb_directory(amoptd['ensembles'], single=False):
        msg = "Cannot import ensembles from the directory: {0}".format(amoptd['ensembles'])
        exit_util.exit_error(msg)

    logger.info("Importing ensembles from directory: {0}".format(amoptd['ensembles']))

    ensembles = glob.glob(os.path.join(amoptd['ensembles'], '*.pdb'))
    amoptd['ensembles'] = ensembles

    # get the data on the ensemble
    ensembles_data = []
    for e in ensembles:
        d = {}
        d['name'] = os.path.splitext(os.path.basename(e))[0]
        d['ensemble_pdb'] = e

        # Get data on the models
        hierarchy = iotbx.pdb.pdb_input(file_name=e).construct_hierarchy()
        d['subcluster_num_models'] = len(hierarchy.models())
        d['num_residues'] = len(hierarchy.models()[0].chains()[0].residue_groups())
        d['ensemble_num_atoms'] = len(hierarchy.models()[0].atoms())

        ensembles_data.append(d)

    amoptd['ensembles_data'] = ensembles_data

    return ensembles


def reorder_models(models, ordered_list_file):
    """Reorder the list of models from a list of models (possibly in a different directory)

    Parameters
    ----------
    models : list, tuple
       A list of model file paths
    ordered_list_file : str
       The path to a file with a model order list

    Returns
    -------
    list
       The model list reordered

    """
    with open(ordered_list_file) as f:
        ordered_list = [line.strip() for line in f if line.strip()]
    mdir = os.path.dirname(models[0])
    model_names = set([os.path.basename(m) for m in models])
    reordered_models = []
    for tm_model in ordered_list:
        name = os.path.basename(tm_model)
        assert name in model_names, "TM model {0} is not in models!".format(name)
        reordered_models.append(os.path.join(mdir, name))
    return reordered_models


def set_phaser_rms_from_subcluster_score(optd):
    """Set the phaser_rms for each ensemble based on its score"""
    if optd['ensemble_options'] is None:
        optd['ensemble_options'] = {}
    for ensemble_data in optd['ensembles_data']:
        name = ensemble_data['name']
        if name not in optd['ensemble_options']:
            optd['ensemble_options'][name] = {}
        optd['ensemble_options'][name]['phaser_rms'] = ensemble_data['subcluster_score']
    return


def sort_ensembles(ensemble_pdbs, ensembles_data=None, keys=None, prioritise=True):
    """Sort AMPLE ensemble data.

    Parameters
    ----------
    ensemble_pdbs: list, tuple
       A list of ensemble file paths
    ensembles_data : list, tuple, optional
       A list of ensembles' data dictionaries
    keys : list, tuple, optional
       A list of keys in ensembles_data by which we sort. The key order must be high > middle > low!
    prioritise : bool, optional
       Favour ensembles in truncation level 20-50% [default: True]

    Returns
    -------
    list
       A sorted list of the ensemble pdb file paths

    Raises
    ------
    RuntimeError
       Not enough ensembles provided

    """
    if len(ensemble_pdbs) < 1:
        raise RuntimeError("Not enough ensembles")
    # Sort with out method only if we have data which contains our generated data
    if ensembles_data:
        ensemble_pdbs_sorted = _sort_ensembles(ensemble_pdbs, ensembles_data, keys, prioritise)
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)
    return ensemble_pdbs_sorted


def _sort_ensembles(ensemble_pdbs, ensemble_data, keys, prioritise):
    """Sort ensembles based on ``keys`` or the following order
        1) Cluster
        2) Truncation level - favoured region of 20-50% first
        3) Subcluster radius threshold
        4) Side chain treatment
    """
    assert len(ensemble_pdbs) == len(ensemble_data), "Unequal ensembles data for sorting"

    # Keys we want to sort data by - differs for different source of ensembles
    if keys is None:
        keys = [
            'cluster_num', 'truncation_score_key', 'truncation_level', 'subcluster_radius_threshold',
            'side_chain_treatment'
        ]
    else:
        keys = [key for key in keys if key in ensemble_data[0] and ensemble_data[0][key]]

    if keys:
        # Zip the data so order remains identical between pdbs and data
        ensembles_zipped = zip(ensemble_pdbs, ensemble_data)

        # The sorting itself
        ensembles_zipped_ordered = sorted(ensembles_zipped, key=lambda e: [e[1][k] for k in keys])

        # Sort the ensembles further if we want to prioritise the "sweet spot"
        # between 20-50% truncation level
        if "truncation_level" in keys and prioritise:
            ensembles_zipped_ordered = _sweet_spotting(ensembles_zipped_ordered, keys)

        # We need to `unzip` the data to get a list of ordered pdbs
        ensemble_pdbs_sorted, _ = zip(*ensembles_zipped_ordered)
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)

    return list(ensemble_pdbs_sorted)


def _sweet_spotting(ensembles_zipped_ordered, keys):
    """Sort the ensembles further if we want to prioritise the "sweet spot"
    between 20-50% truncation level

    Parameters
    ----------
    ensembles_zipped_ordered : list, tuple
       Zipped and __sorted__ data
    keys : list, tuple
       Sorting keys

    Returns
    -------
    list
       Sweet-spotted data

    """

    class Eimer(object):
        """Container unit for binnable data"""

        def __init__(self):
            self.low = []
            self.mid = []
            self.high = []

        def add(self, e):
            if e[1]['truncation_level'] > 50:
                self.high.append(e)
            elif e[1]['truncation_level'] < 20:
                self.low.append(e)
            else:
                self.mid.append(e)

    def magic(key, data):
        """Magic binning happening here"""
        container = {k[key]: [] for k in zip(*data)[1]}
        for e in data:
            container[e[1][key]].append(e)
        mulleimer = []
        for c in sorted(list(container)):
            bin = Eimer()
            for e in container[c]:
                bin.add(e)
            mulleimer.extend(bin.mid + bin.low + bin.high)
        return mulleimer

    # We need to bin the data here in case we want the cluster number or truncation score key to take priority
    # over the truncation level. The binning will, thus, respect the cluster number by binning it accordingly.
    if "cluster_num" in keys and keys.index("cluster_num") < keys.index("truncation_level"):
        return magic("cluster_num", ensembles_zipped_ordered)

    elif "truncation_score_key" in keys and keys.index("truncation_score_key") < keys.index("truncation_level"):
        return magic("truncation_score_key", ensembles_zipped_ordered)

    else:
        bin = Eimer()
        for ens in ensembles_zipped_ordered:
            bin.add(ens)
        return bin.mid + bin.low + bin.high
