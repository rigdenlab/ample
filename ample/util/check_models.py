#!/usr/bin/env ccp4-python
__author__ = "Jens Thomas"

import glob
import logging
import os

import iotbx.pdb

from ample.util import pdb_edit

logger = logging.getLogger(__name__)


class CheckModelsResult():
    def __init__(self):
        self.created_updated_models = False
        self.error = None
        self.homologs = False
        self.nmr = False
        self.merged_chains = False
        self.models_dir = None
        self.single_structure = False

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def check_models_dir(models_in_dir, models_out_dir):
    """Examine a directory of PDB files to determine their suitablilty for running with AMPLE."""
    assert os.path.isdir(models_in_dir)
    pdb_structures = glob.glob(os.path.join(models_in_dir, "*.pdb"))
    pdb_structures += glob.glob(os.path.join(models_in_dir, "*.PDB"))
    results = CheckModelsResult()
    results.models_dir = models_out_dir
    
    results = check_models(pdb_structures, results)
    if not results.created_updated_models:
        # if the models were ok, we can just use as is
        results.models_dir = models_in_dir

    return results
    

def check_models(pdb_structures, results):
    """Examine PDB files to determine their suitablilty for running with AMPLE."""
    assert len(pdb_structures) > 0
    assert results.models_dir is not None
    updated = False
    hierarchies = []
    ref_data = None
    num_pdbs = len(pdb_structures)
    if num_pdbs == 1:
        pdb = pdb_structures[0]
        hierarchy = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        logger.debug("Processing a single input PDB file: {}".format(pdb))
        if len(hierarchy.models()) > 1:
            # Assume NMR model so make sure all have 1 chain with the same sequence
            if not single_chain_models_with_same_sequence(hierarchy):
                results.error = \
                "Supplied with a single pdb file that contained multiple models with multiple or unmatching chains: {}"\
                    .format(pdb)
                return results
            logger.debug("Found a single pdb with multiple models - assuming an NMR ensemble")
            results.nmr = True
        else:
            logger.info("check_models found a single pdb with a single model")
            if not len(hierarchy.only_model().chains()) == 1:
                results.error = \
                    "Supplied with a single pdb file that contained a single model with multiple chains: {}"\
                    .format(pdb)
                return results
            results.single_structure = True
        updated = pdb_edit.add_missing_single_chain_ids(hierarchy)
        if updated:
            hierarchies.append(hierarchy)
    else:
        # multiple pdb structures - get a list of chain_ids and sequences for all members
        multiple = [False] * num_pdbs # tracks if there are multiple chains in the models
        chains_match_across_pdbs = [False] * num_pdbs # do chains with the same index have the same sequence across pdbs?
        logger.debug("Processing {} input PDB files".format(num_pdbs))
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
        # We now have multiple structures, each with one chain - make sure they all have a chain.id
        chains_updated = pdb_edit.add_missing_single_chain_ids(hierarchies)
        updated = chains_updated | updated # see if anything has been updated

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
