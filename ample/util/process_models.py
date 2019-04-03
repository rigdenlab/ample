#!/usr/bin/env ccp4-python
__author__ = "Jens Thomas"

import glob
import logging
import os
import shutil
import tarfile
import zipfile

import iotbx.pdb

from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import exit_util

logger = logging.getLogger(__name__)


class CheckModelsResult():
    def __init__(self):
        self.created_updated_models = False
        self.error = None
        self.homologs = False
        self.ensemble = False
        self.merged_chains = False
        self.models_dir = None
        self.num_structures = 0
        self.num_models = 0
        self.sequence = None
        
    @property
    def single_ensemble(self):
        return self.num_structures == 1 and self.num_models > 1 and self.ensemble

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def extract_and_validate_models(amoptd):
    """Extract models given to AMPLE from arguments in the amoptd and validate
    that they are suitable

    Parameters
    ----------
    amoptd : dict
       AMPLE options dictionary
    """
    
    def is_pdb_file(filename):
        return any(map(lambda x: filename.endswith(x), pdb_suffixes))

    def path_to_quark_alldecoy(pdb_files):
        QUARK_DECOY_NAME = 'alldecoy.pdb'
        for pdb in pdb_files:
            if os.path.basename(pdb) == QUARK_DECOY_NAME:
                return pdb
        return None

    models_arg = amoptd['models'] # command-line arg from user
    if models_arg is None:
        return
    models_dir_final = amoptd['models_dir'] # the directory where the models need to end up
    if models_dir_final is None:
        models_dir_final = os.path.join(amoptd['work_dir'], 'models')
        logger.debug("Setting models_dir_final to: %s", models_dir_final)
        amoptd['models_dir'] = models_dir_final
    pdb_suffixes = ['.pdb', '.PDB']
    quark_models = False
    num_quark_models = 0

    if os.path.isfile(models_arg):
        filepath = models_arg
        models_dir_tmp = os.path.join(amoptd['work_dir'], 'models.tmp')
        if not os.path.isdir(models_dir_tmp):
            os.mkdir(models_dir_tmp)
        # Extract /copy any pdb files into models_dir_tmp
        if tarfile.is_tarfile(filepath):
            pdb_files = ample_util.extract_tar(filepath, models_dir_tmp, suffixes=pdb_suffixes)
        elif zipfile.is_zipfile(filepath):
            pdb_files = ample_util.extract_zip(filepath, models_dir_tmp, suffixes=pdb_suffixes)
        elif is_pdb_file(filepath):
            shutil.copy2(filepath, models_dir_tmp)
            pdb_files = [os.path.join(models_dir_tmp, filepath)]
        else:
            raise RuntimeError("Do not know how to handle input models file: {}".format(filepath))
        # See if we have an alldecoy.pdb file
        quark_decoy = path_to_quark_alldecoy(pdb_files)
        if quark_decoy:
            # Quark decoys are processed by us so go straight into final directory without checking
            num_quark_models = split_quark_alldecoy(quark_decoy, models_dir_final)
            quark_models = True
            amoptd['quark_models'] = quark_models
            shutil.rmtree(models_dir_tmp) # delete as contains uneeded files extracted from archive
    elif os.path.isdir(models_arg):
        models_dir_tmp = models_arg
    elif isinstance(models_arg, list):
        # Assume all models are in the same directory
        models_dir = os.path.dirname(models_arg[0])

    if quark_models:
        # Null result - we extracted the models so assume are ok
        results = CheckModelsResult()
        results.num_structures = num_quark_models
        results.num_models = num_quark_models
        results.models_dir = models_dir_final
    else:
        results = check_models_dir(models_dir_tmp, models_dir_final)

    amoptd['models_dir'] = results.models_dir
    amoptd['models'] = glob.glob(os.path.join(results.models_dir, "*.pdb"))
    return results

def handle_model_import(amoptd, results):
    """Handle any errors flagged up by importing the models."""
    error_msg = None
    if results.error:
        error_msg = "Error importing models: {}".format(results.error)
    elif results.homologs and not amoptd['homologs']:
        error_msg = "Imported models were not sequence identical, but homologs mode wasn't selected"
    if error_msg:
        exit_util.exit_error(error_msg)
    
    if results.single_ensemble and amoptd['webserver_uri']:
        logger.info("** Webserver mode got single NMR model so turning on NMR mode **")
        amoptd['nmr_model_in'] = amoptd['models']


def check_models_dir(models_in_dir, models_out_dir):
    """Examine a directory of PDB files to determine their suitability for running with AMPLE."""
    assert os.path.isdir(models_in_dir)
    pdb_structures = glob.glob(os.path.join(models_in_dir, "*.pdb"))
    pdb_structures += glob.glob(os.path.join(models_in_dir, "*.PDB"))
    results = CheckModelsResult()
    results.models_dir = models_out_dir
    
    results = check_models(pdb_structures, results)
    if not results.created_updated_models:
        # if the models were ok, we can just use as is
        results.models_dir = models_in_dir
    logger.info("Using models from directory: %s" % results.models_dir )
    return results
    

def check_models(pdb_structures, results):
    """Examine PDB files to determine their suitablilty for running with AMPLE."""
    assert len(pdb_structures) > 0
    assert results.models_dir is not None
    updated = False
    hierarchies = []
    ref_data = None
    num_structures = len(pdb_structures)
    if num_structures == 1:
        pdb = pdb_structures[0]
        logger.debug("Processing a single input PDB file: {}".format(pdb))
        hierarchy = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        num_models = len(hierarchy.models())
        results.num_models = num_models
        if num_models > 1:
            if not single_chain_models_with_same_sequence(hierarchy):
                results.error = \
                "Supplied with a single pdb file that contained multiple models with multiple or unmatching chains: {}"\
                    .format(pdb)
                return results
            logger.debug("Found a single pdb with multiple models - assuming an NMR ensemble")
            results.ensemble = True
            results.sequence = pdb_edit.chain_sequence(hierarchy.models()[0].only_chain())
        else:
            logger.info("check_models found a single pdb with a single model")
            if not len(hierarchy.only_model().chains()) == 1:
                results.error = \
                    "Supplied with a single pdb file that contained a single model with multiple chains: {}"\
                    .format(pdb)
                return results
        updated = pdb_edit.add_missing_single_chain_ids(hierarchy)
        if updated:
            hierarchies.append(hierarchy)
    else:
        # multiple pdb structures - get a list of chain_ids and sequences for all members
        multiple = [False] * num_structures # tracks if there are multiple chains in the models
        chains_match_across_pdbs = [False] * num_structures # do chains with the same index have the same sequence across pdbs?
        logger.debug("Processing {} input PDB files".format(num_structures))
        for idx_pdb, pdb in enumerate(pdb_structures):
            h = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
            if len(h.models()) > 1:
                # Assume an NMR ensemble so just extract the first member as representative of all
                logger.info("Multiple models in pdb {}. Assuming NMR ensemble so extracting first model".format(pdb))
                h = iotbx.pdb.hierarchy.root()
                h.append_model((h.models()[0].detached_copy()))
                updated = True
            hierarchies.append(h)
            seq_data = pdb_edit._sequence_data(h)
            if ref_data is None:
                # Get sequence data for first pdb
                ref_data = seq_data
                if len(ref_data) > 1:
                    multiple[idx_pdb] = True
                continue
            else:
                check_sequences_match(ref_data, seq_data, idx_pdb, chains_match_across_pdbs, multiple)
            # We now have a list of single-model structures
            results.num_models = num_structures
        if any(multiple):
            logger.debug("Processing multichain pdbs")
            if sum(multiple) != len(multiple):
                results.error = "check_models: given multiple models, but not all had multiple chains"
                return results
            if sum(chains_match_across_pdbs) != len(chains_match_across_pdbs):
                results.error = "check_models: given multiple models with multiple chains, but chains did not match across pdbs"
                return results
            # merge all chains
            for i, h in enumerate(hierarchies):
                hierarchies[i] = pdb_edit._merge_chains(h)
            results.merged_chains = True
            updated = True
        else:
            # multiple single-chain pdbs - are they homologs?
            logger.debug("Processing single-chain pdbs")
            if sum(chains_match_across_pdbs) != len(chains_match_across_pdbs):
                logger.debug("Chains don't match across pdbs so assuming homologs")
                results.homologs = True
            else:
                # All have the same sequence
                results.sequence = pdb_edit.chain_sequence(hierarchies[0].only_model().only_chain())
        # We now have multiple structures, each with one chain - make sure they all have a chain.id
        chains_updated = pdb_edit.add_missing_single_chain_ids(hierarchies)
        updated = chains_updated | updated # see if anything has been updated

    results.num_structures = num_structures
    if updated:
        # write out new files
        results.created_updated_models = True
        for pdb, h in zip(pdb_structures, hierarchies):
            basename = os.path.basename(pdb)
            pdbout = os.path.join(results.models_dir, basename)
            with open(pdbout, 'w') as f:
                f.write("REMARK Original file:{}\n".format(pdb))
                f.write(h.as_pdb_string(anisou=False))  
    return results


def check_sequences_match(ref_data, seq_data, idx_pdb, chains_match_across_pdbs, multiple):
    """Check whether chains with the same index across multiple pdbs have matching sequences."""
    all_chains_match = 0
    for i, (rd, sd) in enumerate(zip(ref_data, seq_data)):
        ref_seq = ref_data[rd][0]
        seq = seq_data[sd][0]
        if ref_seq == seq:
            all_chains_match += 1
    if all_chains_match == len(ref_data):
        chains_match_across_pdbs[idx_pdb] = True
        if idx_pdb == 1:
            # if the second chain matches the first then the first matches 
            chains_match_across_pdbs[0] = True
    if i > 0:
        multiple[idx_pdb] = True


def single_chain_models_with_same_sequence(hierarchy):
    """Return True if hierarchy contains sequence-identical single-chain models."""
    root_seq = None
    for model in hierarchy.models():
        try:
            chain = model.only_chain()
        except AssertionError:
            return False
        seq = pdb_edit.chain_sequence(chain)
        if root_seq is None:
            root_seq = seq
            continue
        if seq != root_seq:
            return False
    return True

def split_quark_alldecoy(alldecoy, directory):
    """Split a single QUARK PDB with multiple models into individual PDB files

    Parameters
    ----------
    alldecoy : str
       Single QUARK PDB file with multiple model entries
    directory : str
       Directory to extract the PDB files to

    Returns
    -------
    extracted_models : list
       List of PDB files for all models

    """
    logger.info("Extracting decoys from: %s into %s", alldecoy, directory)
    if not os.path.isdir(directory):
        os.mkdir(directory)
    smodels = []
    with open(alldecoy, 'r') as f:
        m = []
        for line in f:
            if line.startswith("ENDMDL"):
                m.append(line)
                smodels.append(m)
                m = []
            else:
                m.append(line)

    if not len(smodels):
        raise RuntimeError("Could not extract any models from: {0}".format(alldecoy))

    for i, m in enumerate(smodels):
        fpath = os.path.join(directory, "quark_{0}.pdb".format(i))
        with open(fpath, 'w') as f:
            for line in m:
                #  Reconstruct something sensible as from the coordinates on it's all quark-specific
                # and there is no chain ID
                if line.startswith("ATOM"):
                    line = line[:21] + 'A' + line[22:54] + "  1.00  0.00              \n"
                f.write(line)
            logger.debug("Wrote: %s", fpath)
    num_models = i + 1
    logger.info("Extracted %d models from quark archive" % num_models)
    return num_models
