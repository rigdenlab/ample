"""Ensembler module for single model structures"""

__author__ = "Felix Simkovic, Jens Thomas and Adam Simpkin"
__date__ = "16 Feb 2016"
__version__ = "1.0"

import gemmi
import logging
import math
import os
import pandas as pd
import string
import sys

from ample.ensembler import _ensembler, truncation_util
from ample.ensembler.constants import SIDE_CHAIN_TREATMENTS
from ample.util import ample_util, pdb_edit

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

    def generate_ensembles(
        self,
        models,
        ensembles_directory=None,
        nproc=None,
        percent_truncation=None,
        percent_fixed_intervals=None,
        side_chain_treatments=SIDE_CHAIN_TREATMENTS,
        truncation_method=None,
        truncation_pruning=None,
        truncation_scorefile=None,
        truncation_scorefile_header=None,
    ):
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

        # standardise the structure
        std_models_dir = os.path.join(self.work_dir, "std_models")
        os.mkdir(std_models_dir)

        std_model = ample_util.filename_append(models[0], 'std', std_models_dir)

        #Clean up potential oddities in input files i.e. missing chain ID, odd PDB headers
        struct = gemmi.read_structure(models[0])
        model = struct[0]
        alphabet = list(string.ascii_lowercase)
        for i, chain in enumerate(model):
            if chain.name == "":
                chain.name = alphabet[i]
        for chain in model:
            for idx, residue in enumerate(chain):
                residue.seqid.num = idx + 1
        pdb_string = [line for line in struct.make_minimal_pdb().split('\n')]
        with open("tmp.pdb", "w") as f_out:
            for line in pdb_string:
                f_out.write(line + os.linesep)

        pdb_edit.standardise(pdbin="tmp.pdb", pdbout=std_model)
        os.unlink("tmp.pdb")

        if truncation_method == truncation_util.TRUNCATION_METHODS.ERRORS:
            self._modify_bfactors(std_model, std_model)

        std_models = [std_model]
        logger.info('Standardised input model: %s', std_models[0])

        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory):
            os.mkdir(self.ensembles_directory)

        truncate_dir = os.path.join(self.work_dir, "single_truncate")
        if not os.path.isdir(truncate_dir):
            os.mkdir(truncate_dir)

        if truncation_method == truncation_util.TRUNCATION_METHODS.SCORES:
            if len(truncation_scorefile_header) < 2:
                msg = "At least two header options for scorefile are required"
                logger.critical(msg)
                raise RuntimeError(msg)

            # Read all the scores into a per residue dictionary
            assert len(truncation_scorefile_header) > 1, "At least two column labels are required"
            residue_scores = self._read_scorefile(truncation_scorefile)
            residue_key = truncation_scorefile_header.pop(0)
            truncation_scorefile_header = list(map(str.strip, truncation_scorefile_header))
            assert all(
                h in residue_scores[0] for h in truncation_scorefile_header
            ), "Not all column labels are in your CSV file"
            self.ensembles = []
            for score_key in truncation_scorefile_header:
                self._generate_ensembles(residue_key, score_key, residue_scores, truncate_dir, std_models,
                                         truncation_method, percent_truncation, percent_fixed_intervals,
                                         truncation_pruning, side_chain_treatments)

        if truncation_method in [truncation_util.TRUNCATION_METHODS.BFACTORS, truncation_util.TRUNCATION_METHODS.ERRORS]:
            struct = gemmi.read_structure(std_model)
            residue_scores = []
            residue_key = "Residue"
            score_key = "Bfactor"
            for chain in struct[0]:
                for residue in chain:
                    residue_scores.append({residue_key: residue.seqid.num, score_key: residue[0].b_iso})
            self._generate_ensembles(residue_key, score_key, residue_scores, truncate_dir, std_models,
                                     truncation_method, percent_truncation, percent_fixed_intervals,
                                     truncation_pruning, side_chain_treatments)

        return self.ensembles

    def _generate_ensembles(self, residue_key, score_key, residue_scores, truncate_dir, std_models, truncation_method,
                            percent_truncation, percent_fixed_intervals, truncation_pruning, side_chain_treatments):
        zipped_scores = self._generate_residue_scorelist(residue_key, score_key, residue_scores)
        score_truncate_dir = os.path.join(truncate_dir, "{}".format(score_key))
        if not os.path.isdir(score_truncate_dir):
            os.mkdir(score_truncate_dir)
        self.ensembles = []
        self.truncator = truncation_util.Truncator(work_dir=score_truncate_dir)
        self.truncator.theseus_exe = self.theseus_exe
        for truncation in self.truncator.truncate_models(
                models=std_models,
                truncation_method=truncation_method,
                percent_truncation=percent_truncation,
                percent_fixed_intervals=percent_fixed_intervals,
                truncation_pruning=truncation_pruning,
                residue_scores=zipped_scores,
        ):

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

            for ensemble in self.edit_side_chains(pre_ensemble, side_chain_treatments, single_structure=True):
                self.ensembles.append(ensemble)

    def generate_ensembles_from_amoptd(self, models, amoptd):
        """Generate ensembles from data in supplied ample data dictionary."""
        kwargs = {
            'percent_truncation': amoptd['percent'],
            'percent_fixed_intervals': amoptd['percent_fixed_intervals'],
            'side_chain_treatments': amoptd['side_chain_treatments'],
            'truncation_method': amoptd['truncation_method'],
            'truncation_pruning': amoptd['truncation_pruning'],
            'truncation_scorefile': amoptd['truncation_scorefile'],
            'truncation_scorefile_header': amoptd['truncation_scorefile_header'],
        }
        if sys.version_info.major == 3:
            kwargs = {k: v for k, v in kwargs.items() if v is not None}
        else:
            kwargs = {k: v for k, v in kwargs.iteritems() if v is not None}
        return self.generate_ensembles(models, **kwargs)

    @staticmethod
    def _modify_bfactors(pdbin, pdbout):
        """Modify error estimates to B-factors"""
        multiplier = 8.0 / 3.0 * math.pi ** 2
        bmax = 999.0
        rms_max = math.sqrt(bmax / multiplier)
        rms_big = 4.0

        struct = gemmi.read_structure(pdbin)
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        occ = max(min(1.0 - (atom.b_iso - rms_big) / (rms_max - rms_big), 1.0), 0.0)
                        bfac = min(multiplier * atom.b_iso ** 2, bmax)
                        atom.occ = occ
                        atom.b_iso = bfac

        pdb_string = [line for line in struct.make_minimal_pdb().split('\n') if not line.startswith('ANISOU')]
        with open(pdbout, "w") as f_out:
            for line in pdb_string:
                f_out.write(line + os.linesep)

    @staticmethod
    def _generate_residue_scorelist(residue_key, score_key, scores):
        """Generate a zipped list of residue indexes and corresponding scores
        
        :residue_key: residue column header keyword
        :score_key: score column header keyword
        :scores: list of dictionaries for each residue
        
        :returns: zipped list of residue index plus score
        """
        assert residue_key in scores[0], "Cannot find residue key {} in scoresfile header: {}".format(
            residue_key, scores[0]
        )
        assert score_key in scores[0], "Cannot find score key {} in scoresfile header: {}".format(score_key, scores[0])
        return [(i[residue_key], i[score_key]) for i in scores]

    @staticmethod
    def _read_scorefile(scorefile):
        """
        :scorefile: CSV score file INCLUDING header line
        
        :returns: list of per residue dictionaries containing column data
        """
        df = pd.read_csv(scorefile)
        df.rename(columns=lambda x: x.strip(), inplace=True)
        return list(df.T.to_dict().values())
