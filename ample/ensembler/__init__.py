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
import os
import shutil
import sys

from constants import SIDE_CHAIN_TREATMENTS
from abinitio import AbinitioEnsembler
from homologs import HomologEnsembler
from single_model import SingleModelEnsembler
from ample.util import ample_util
from ample.util import argparse_util
from ample.util import exit_util
from ample.util import pdb_edit
from ample.util import printTable

logger = logging.getLogger(__name__)

def add_argparse_options(parser, standalone=False):
    """Function to add the ensemble-specific options
    
    Parameters
    ----------
    parser : :obj:`argparse.ArgumentParser`
    standalone : bool, optional
       Define if the parser is standalone [default: False]

    """
    
    if standalone:
        # Add options for running as a standalone module
        argparse_util.add_general_options(parser)
    else:
        parser = parser.add_argument_group('Ensemble Options')

    # Ensenble-specific options here
    parser.add_argument('-cluster_dir',
                       help='Path to directory of pre-clustered models to import')
    
    parser.add_argument('-cluster_method',
                       help='How to cluster the models for ensembling (spicker|fast_protein_cluster')

    parser.add_argument('-ensembler_timeout', type=int,
                       help='Time in seconds before timing out ensembling')

    parser.add_argument('-gesamt_exe', metavar='gesamt_exe',
                       help='Path to the gesamt executable')

    parser.add_argument('-homologs', metavar='True/False',
                       help='Generate ensembles from homologs models (requires -alignment_file)')
    
    parser.add_argument('-homolog_aligner', metavar='homolog_aligner',
                       help='Program to use for structural alignment of homologs (gesamt|mustang)')

    parser.add_argument('-ensemble_max_models',
                       help='Maximum number of models permitted in an ensemble')

    parser.add_argument('-maxcluster_exe',
                       help='Path to Maxcluster executable')
    
    parser.add_argument('-mustang_exe', metavar='mustang_exe',
                       help='Path to the mustang executable')
    
    parser.add_argument('-num_clusters', type=int,
                       help='The number of Spicker clusters of the original decoys that will be sampled [1]')

    parser.add_argument('-percent', metavar='percent_truncation',
                       help='percent interval for truncation')
            
    parser.add_argument('-score_matrix',
                       help='Path to score matrix for spicker')
    
    parser.add_argument('-score_matrix_file_list',
                       help='File with list of ordered model names for the score_matrix')
     
    parser.add_argument('-side_chain_treatments', type=str, nargs='+', action='append',
                        help='The side chain treatments to use. Default: {0}'.format(SIDE_CHAIN_TREATMENTS))

    parser.add_argument('-subcluster_radius_thresholds', type=float, nargs='+',
                        help='The radii to use for subclustering the truncated ensembles')
    
    parser.add_argument('-subcluster_program',
                        help='Program for subclustering models [maxcluster]')

    parser.add_argument('-theseus_exe', metavar='Theseus exe (required)',
                        help='Path to theseus executable')
    
    parser.add_argument('-thin_clusters', metavar='True/False',
                        help='Create ensembles from 10 clusters with 1 + 3A subclustering and polyAlanine sidechains')

    parser.add_argument('-truncation_method',
                        help='How to truncate the models for ensembling percent|thresh|focussed|scores')

    parser.add_argument('-truncation_pruning',
                        help='Whether to remove isolated residues (single)')

    parser.add_argument('-truncation_scorefile',
                        help="CSV file containing per residue scores - COLUMN ONE MUST BE RESIDUE INDEX STARTING FROM 1")

    parser.add_argument('-truncation_scorefile_header', nargs='+',
                        help="column headers to be used to create ensembles")

    return
                                                                                                                                                                   

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
    if not (amoptd['single_model'] or amoptd['models_dir']):
        msg = 'AMPLE ensembler needs either a single_model or a models_dir argument'
        exit_util.exit_error(msg, sys.exc_info()[2])
        if amoptd['single_model'] and not os.path.isfile(amoptd['single_model']):
            msg = 'Cannot find single_model pdb: {0}'.format(amoptd['single_model'])
            exit_util.exit_error(msg, sys.exc_info()[2])
        elif amoptd['models_dir'] and not os.path.isfile(amoptd['models_dir']):
            msg = 'Cannot find models_dir: {0}'.format(amoptd['models_dir'])
            exit_util.exit_error(msg, sys.exc_info()[2])
            
    models = list([amoptd['single_model']]) if amoptd['single_model_mode'] else \
        glob.glob(os.path.join(amoptd['models_dir'], "*.pdb"))

    #if amoptd['cluster_method'] == 'spicker_tmscore':
    #    models = reorder_models(models, amoptd['score_matrix_file_list'])
    #    if not (os.path.isfile(amoptd['score_matrix']) and os.path.isfile(amoptd['score_matrix_file_list'])):
    #        raise RuntimeError("spicker_tmscore needs a score_matrix and score_matrix_file_list")
    #    ensembler.score_matrix = amoptd['score_matrix']

    # Run ensemble creation
    ensembles = ensembler.generate_ensembles_from_amoptd(models, amoptd)
    
    ############################################################################
    # Hack to pull out the data - need to update code to work with ensemble objects rather than dictionaries
    amoptd['ensembles'] = [e.pdb for e in ensembles]
    amoptd['ensembles_data'] = [e.__dict__ for e in ensembles]

    # We need to let the main process know that we have succeeded as this module could be run on a cluster node with no link
    # to the parent process, so we create a file here indicating that we got this far and didn't die from an exception
    with open(amoptd['ensemble_ok'],'w') as f: f.write('ok\n')

    # Delete all intermediate files if we're purging
    if amoptd['purge']: shutil.rmtree(ensembler.work_dir)
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
    clusters : dict
    cluster_method : str
    truncation_method : str
    truncation_percent : float
    side_chain_treatments : list
    
    """
    clusters = {}  # Loop through all ensemble data objects and build up a data tree
    cluster_method = None
    truncation_method = None
    truncation_percent = None
    side_chain_treatments = []
    for e in ensembles_data:
        if not cluster_method:
            cluster_method = e['cluster_method']
            truncation_percent = e['truncation_percent']
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
        if sct not in side_chain_treatments: side_chain_treatments.append(sct)
        if sct not in clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct']:
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct] = {}
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['name'] = e['name']
            clusters[cnum]['tlevels'][tlvl]['radius_thresholds'][srt]['sct'][sct]['num_atoms'] = e['ensemble_num_atoms']

    return clusters, cluster_method, truncation_method, truncation_percent, side_chain_treatments


def cluster_table_data(clusters, cluster_num, side_chain_treatments):
    """Generate a table containing the cluster data
    
    Parameters
    ----------
    clusters : list, tuple 
       A list of cluster dictionaries
    cluster_num : int
       The cluster number of interest
    side_chain_treatments : list, tuple
       A list of side chain treatments

    Returns
    -------
    list
       A list of table data lines

    """
    # FIX TO IGNORE side_chain_treatments
    # tdata = [("Name", "Truncation Level", u"Variance Threshold (\u212B^2)", "No. Residues", u"Radius Threshold (\u212B)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
    tdata = [("Name", "Cluster", "Truncation Level", "Variance Threshold (A^2)", "No. Residues", "Radius Threshold (A)", "No. Decoys", "Number of Atoms", "Sidechain Treatment")]
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
    clusters, cluster_method, truncation_method, truncation_percent, side_chain_treatments = collate_cluster_data(ensembles_data)
    num_clusters = len(clusters)

    tableFormat = printTable.Table()
    rstr = "\n"
    rstr += "Ensemble Results\n"
    rstr += "----------------\n\n"
    rstr += "Cluster method: {0}\n".format(cluster_method)
    rstr += "Truncation method: {0}\n".format(truncation_method)
    rstr += "Percent truncation: {0}\n".format(truncation_percent   )
    rstr += "Number of clusters: {0}\n".format(num_clusters)

    for cluster_num in sorted(clusters.keys()):
        rstr += "\n"
        rstr += "Cluster {0}\n".format(cluster_num)
        rstr += "Number of models: {0}\n".format(clusters[cluster_num]['cluster_num_models'])
        rstr += "Cluster centroid: {0}\n".format(clusters[cluster_num]['cluster_centroid'])
        rstr += "\n"
        tdata = cluster_table_data(clusters, cluster_num, side_chain_treatments)
        rstr += tableFormat.pprint_table(tdata)

    rstr += "\nGenerated {0} ensembles\n\n".format(len(ensembles_data))

    return rstr


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
        ordered_list = [ line.strip() for line in f if line.strip() ]
    mdir = os.path.dirname(models[0])
    model_names = set([ os.path.basename(m) for m in models ])
    reordered_models = []
    for tm_model in ordered_list:
        name = os.path.basename(tm_model)
        assert name in model_names,"TM model {0} is not in models!".format(name)
        reordered_models.append(os.path.join(mdir, name))
    return reordered_models


def set_phaser_rms_from_subcluster_score(optd):
    """Set the phaser_rms for each ensemble based on its score"""
    if optd['ensemble_options'] is None: optd['ensemble_options'] = {}
    for ensemble_data in optd['ensembles_data']:
        name = ensemble_data['name']
        if name not in optd['ensemble_options']: optd['ensemble_options'][name] = {}
        optd['ensemble_options'][name]['phaser_rms'] = ensemble_data['subcluster_score']
    return


def sort_ensembles(ensemble_pdbs, ensembles_data=None, prioritise=True):
    """Sort AMPLE ensemble data.
    
    Parameters
    ----------
    ensemble_pdbs: list 
       A list of ensemble file paths
    ensembles_data : list, tuple, optional
       A list of ensembles' data dictionaries
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
        ensemble_pdbs_sorted = _sort_ensembles(ensemble_pdbs, ensembles_data, prioritise)    
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)
    return ensemble_pdbs_sorted


def _sort_ensembles(ensemble_pdbs, ensembles_data, prioritise):
    """Sort ensembles based on
        1) Cluster
        2) Truncation level - favoured region of 20-50% first
        3) Subcluster radius threshold
        4) Side chain treatment
    """
    assert len(ensemble_pdbs) == len(ensembles_data), "Unequal ensembles data for sorting"

    # Determine which keys we can sort our data with
    def _get_keys(ensemble):
        _criteria = ['cluster_num',
                     'truncation_score_key',
                     'truncation_level',
                     'subcluster_radius_threshold',
                     'side_chain_treatment']

        for crit in _criteria:
            if ensemble.has_key(crit) and ensemble[crit]:
                yield crit

    # Keys we want to sort data by - differs for different source of ensembles
    keys_to_sort = list(_get_keys(ensembles_data[0]))
    if keys_to_sort:
        # Zip the data so order remains identical between pdbs and data
        ensembles_zipped = zip(ensemble_pdbs, ensembles_data)
        if "truncation_level" in keys_to_sort and prioritise:
            ensembles_zipped_ordered = _sort_ensembles_prioritise(ensembles_zipped, keys_to_sort)
        else:
            ensembles_zipped_ordered = _sort_ensembles_parameters(ensembles_zipped, keys_to_sort)
        # We need to `unzip` the data to get a list of ordered pdbs
        ensemble_pdbs_sorted, _ = zip(*ensembles_zipped_ordered)
    else:
        ensemble_pdbs_sorted = sorted(ensemble_pdbs)

    return ensemble_pdbs_sorted


def _sort_ensembles_prioritise(ensembles_zipped, keys_to_sort):
    """Sort our ensembles based on the data in order of the keys provided"""
    # Group the ensembles either by cluster ... or by truncation score ...
    # ... otherwise throw them all in one pot
    tmp_data = collections.defaultdict(list)
    if "cluster_num" in keys_to_sort:
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['cluster_num']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys(), key=int)

    elif "truncation_score_key" in keys_to_sort:
        for ensemble in ensembles_zipped:
            tmp_data[ensemble[1]['truncation_score_key']].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())

    else:
        for ensemble in ensembles_zipped:
            tmp_data["all"].append(ensemble)
        iterator_keys = sorted(tmp_data.keys())

    # Iterate through clusters and group based on truncation level
    ensembles_zipped_ordered = []
    for sbin in iterator_keys:
        low, mid, high = [], [], []
        for ensemble in _sort_ensembles_parameters(tmp_data[sbin], keys_to_sort):
            if ensemble[1]['truncation_level'] > 50:
                high.append(ensemble)
            elif ensemble[1]['truncation_level'] < 20:
                low.append(ensemble)
            else:
                mid.append(ensemble)
        ensembles_zipped_ordered += mid + low + high

    return ensembles_zipped_ordered


def _sort_ensembles_parameters(ensembles_zipped, keys_to_sort):
    """Tiny wrapper function to sort ensembles data based on keys provided
    """
    def _extract(ens):
        return [ens[1][crit] for crit in keys_to_sort]
    return sorted(ensembles_zipped, key=_extract)

