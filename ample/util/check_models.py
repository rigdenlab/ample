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
        self.homolog = False
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
    assert results.models_dir is not None
    updated = False
    hierarchies = []
    if len(pdb_structures) == 1:
        pdb = pdb_structures[0]
        hierarchy = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        if len(hierarchy.models()) > 1:
            # Assume NMR model so make sure all have 1 chain with the same sequence
            if not single_chain_models_with_same_sequence(hierarchy):
                results.error = \
                "Supplied with a single pdb file that contained multiple models with multiple or unmatching chains: {}"\
                    .format(pdb)
                return results
            logger.info("check_models found a single pdb with multiple models - assuming an NMR ensemble")
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
        ref_data = None
        multiple = [False] * len(pdb_structures) # tracks if there are multiple chains in the models
        for idx_pdb, pdb in enumerate(pdb_structures):
            h = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
            if len(h.models) > 1:
                # Assume an NMR ensemble so just extract the first member
                logger.debug("Multiple models in pdb {} so extracting first model".format(pdb))
                h = iotbx.pdb.hierarchy.root()
                h.append_model((h.models()[0].detached_copy()))
                updated = True
            hierarchies.append(h)
            seq_data = pdb_edit._sequence_data(h)
            if ref_data is None:
                ref_data = seq_data
                continue
            for i, (rd, sd) in enumerate(zip(ref_data, seq_data)):
                ref_seq = ref_data[rd][0]
                seq = seq_data[sd][0]
                if ref_seq != seq:
                    results.error = "check_models: pdb {} chain number {} has different sequence from chain {} in first pdb {}"\
                        .format(pdb, i, i, pdb_structures[0])
                    return results
            if i > 0:
                multiple[idx_pdb] = True
        # We now have a list of pdbs with 1 or more chains, but with chain sequences that match across all structures
        if any(multiple):
            if sum(multiple) != len(multiple):
                results.error = "check_models: given multiple models, but not all had multiple chains"
                return results
            # merge all chains
            for i, h in enumerate(hierarchies):
                hierarchies[i] = pdb_edit._merge_chains(h)
            results.merged_chains = True
            updated = True
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
