"""Ensembler module for single model structures"""

__author__ = "Felix Simkovic, and Jens Thomas"
__date__ = "16 Feb 2016"
__version__ = "1.0"

import csv
import logging
import os
import pandas as pd

import _ensembler
import truncation_util
from constants import SIDE_CHAIN_TREATMENTS
from ample.util import ample_util
from ample.util import pdb_edit

logger = logging.getLogger(__name__)


class SingleModelEnsembler(_ensembler.Ensembler):
    """Ensemble creator using on a single input structure and a corresponding
       score file with per residue scores for truncation
    """

    def __init__(self, **kwargs):
        # Inherit all functions from Parent Ensembler
        super(SingleModelEnsembler, self).__init__(**kwargs)

        # Set SingleModelEnsembler specific parameters
        self.truncation_scorefile = None

        return

    def generate_ensembles(self,
                           models,
                           ensembles_directory=None,
                           nproc=None,
                           percent_truncation=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS,
                           truncation_method=None,
                           truncation_pruning=None,
                           truncation_scorefile=None,
                           truncation_scorefile_header=None):
        """Method to generate ensembles from a single structure based on 
        residue scores"""

        if not truncation_method:
            truncation_method = self.truncation_method
        if not truncation_pruning:
            truncation_pruning = self.truncation_pruning
        if not truncation_scorefile:
            truncation_scorefile = self.truncation_scorefile

        if len(models) > 1:
            msg = "More than 1 structure provided"
            logger.critical(msg)
            raise RuntimeError(msg)

        if len(truncation_scorefile_header) < 2:
            msg = "At least two header options for scorefile are required"
            logger.critical(msg)
            raise RuntimeError(msg)

        # standardise the structure
        std_models_dir = os.path.join(self.work_dir, "std_models")
        os.mkdir(std_models_dir)

        std_model = ample_util.filename_append(
            models[0], 'std', std_models_dir)
        pdb_edit.standardise(
            pdbin=models[0], pdbout=std_model, del_hetatm=True)
        std_models = [std_model]
        logger.info('Standardised input model: %s', std_models[0])

        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory):
            os.mkdir(self.ensembles_directory)

        truncate_dir = os.path.join(self.work_dir, "single_truncate")
        if not os.path.isdir(truncate_dir):
            os.mkdir(truncate_dir)

        # Read all the scores into a per residue dictionary
        assert len(
            truncation_scorefile_header) > 1, "At least two column labels are required"
        residue_scores = self._read_scorefile(truncation_scorefile)
        residue_key = truncation_scorefile_header.pop(0).lower()
        assert all(h in residue_scores[0] for h in truncation_scorefile_header), \
            "Not all column labels are in your CSV file"

        self.ensembles = []
        for score_key in truncation_scorefile_header:
            zipped_scores = self._generate_residue_scorelist(residue_key,
                                                             score_key,
                                                             residue_scores)

            score_truncate_dir = os.path.join(
                truncate_dir, "{0}".format(score_key))
            if not os.path.isdir(score_truncate_dir):
                os.mkdir(score_truncate_dir)

            self.truncator = truncation_util.Truncator(
                work_dir=score_truncate_dir)
            self.truncator.theseus_exe = self.theseus_exe
            for truncation in self.truncator.truncate_models(models=std_models,
                                                             truncation_method=truncation_method,
                                                             percent_truncation=percent_truncation,
                                                             truncation_pruning=truncation_pruning,
                                                             residue_scores=zipped_scores):

                # Create Ensemble object
                pre_ensemble = _ensembler.Ensemble()
                pre_ensemble.num_residues = truncation.num_residues
                pre_ensemble.truncation_dir = truncation.directory
                pre_ensemble.truncation_level = truncation.level
                pre_ensemble.truncation_method = truncation.method
                pre_ensemble.truncation_percent = truncation.percent
                pre_ensemble.truncation_residues = truncation.residues
                pre_ensemble.truncation_variance = truncation.variances
                pre_ensemble.truncation_score_key = score_key.lower()
                pre_ensemble.pdb = truncation.models[0]

                for ensemble in self.edit_side_chains(pre_ensemble,
                                                      side_chain_treatments,
                                                      single_structure=True):
                    self.ensembles.append(ensemble)

        return self.ensembles

    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary."""
        kwargs = {'percent_truncation': amoptd['percent'],
                  'side_chain_treatments': amoptd['side_chain_treatments'],
                  'truncation_method': amoptd['truncation_method'],
                  'truncation_pruning': amoptd['truncation_pruning'],
                  'truncation_scorefile': amoptd['truncation_scorefile'],
                  'truncation_scorefile_header': amoptd['truncation_scorefile_header']}
        # strip out any that are None
        kwargs = {k: v for k, v in kwargs.iteritems() if v is not None}
        return self.generate_ensembles(models, **kwargs)

    # staticmethod so that we can test without instantiating an Ensembler
    @staticmethod
    def _generate_residue_scorelist(residue_key, score_key, scores):
        """Generate a zipped list of residue indexes and corresponding scores
        
        :residue_key: residue column header keyword
        :score_key: score column header keyword
        :scores: list of dictionaries for each residue
        
        :returns: zipped list of residue index plus score
        """
        assert residue_key in scores[0], "Cannot find residue key in scoresfile"
        assert score_key in scores[0], "Cannot find score key in scoresfile"
        return [(i[residue_key], i[score_key]) for i in scores]

    # staticmethod so that we can test without instantiating an Ensembler
    @staticmethod
    def _read_scorefile(scorefile):
        """
        :scorefile: CSV score file INCLUDING header line
        
        :returns: list of per residue dictionaries containing column data
        """
        return pd.read_csv(scorefile).T.to_dict().values()
