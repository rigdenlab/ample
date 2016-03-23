#!/usr/bin/env ccp4-python

"""
16.02.2016

@author: hlfsimko
"""

# System
import copy
import csv
import logging
import os

# Custom
from ample.ensembler import _ensembler
from ample.ensembler.constants import SIDE_CHAIN_TREATMENTS
from ample.util import ample_util
from ample.util import pdb_edit

_logger = logging.getLogger(__name__)


class Ensembler(_ensembler.Ensembler):
    """Ensemble creator using on a single input structure and a corresponding
       score file with per residue scores for truncation
    """
        
    def __init__(self):
        # Inherit all functions from Parent Ensembler
        _ensembler.Ensembler.__init__(self)
        
        # Set Ensembler_Single specific parameters
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
                           truncation_scorefile_header=None,
                           work_dir=None):
        """Method to generate ensembles from a single structure based on 
        residue scores"""
        
        # Work dir set each time
        if not work_dir: raise RuntimeError, "Need to set work_dir!"
        self.work_dir = work_dir
        
        if not truncation_method:
            truncation_method = self.truncation_method
        if not truncation_pruning:
            truncation_pruning = self.truncation_pruning
        if not truncation_scorefile:
            truncation_scorefile = self.truncation_scorefile
        if not ensembles_directory:
            self.ensembles_directory = os.path.join(work_dir, "ensembles")
        else:
            self.ensembles_directory = ensembles_directory

        if len(models) > 1:
            msg = "More than 1 structure provided"
            _logger.critical(msg)
            raise RuntimeError(msg)
        
        if len(truncation_scorefile_header) < 2:
            msg = "At least two header options for scorefile are required"
            _logger.critical(msg)
            raise RuntimeError(msg)

        # standardise the structure
        std_models_dir = os.path.join(work_dir, "std_models")
        os.mkdir(std_models_dir)
        
        std_model = ample_util.filename_append(models[0], 'std', std_models_dir)
        pdb_edit.standardise(pdbin=models[0], pdbout=std_model, del_hetatm=True)
        std_models = [std_model]
        _logger.info('Standardised input model: {0}'.format(std_models[0]))
              
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)

        truncate_dir = os.path.join(self.work_dir,"single_truncate")
        if not os.path.isdir(truncate_dir): os.mkdir(truncate_dir)
        
        # Read all the scores into a per residue dictionary
        assert len(truncation_scorefile_header) > 1, "At least two column labels are required"
        residue_scores = self._read_scorefile(truncation_scorefile)
        residue_key = truncation_scorefile_header.pop(0).lower()
        assert all(i.lower() in residue_scores[0].keys() for i in truncation_scorefile_header), \
                "Not all column labels are in your CSV file"
        
        self.ensembles = []
        self.ensembles_data = []
        
        for score_key in truncation_scorefile_header:
            score_key = score_key.lower()
            zipped_scores = self._generate_residue_scorelist(residue_key, 
                                                             score_key, 
                                                             residue_scores)

            score_truncate_dir = os.path.join(truncate_dir, "{0}".format(score_key))
            if not os.path.isdir(score_truncate_dir): os.mkdir(score_truncate_dir)
            
            for truncated_model, truncated_model_data, truncated_model_dir in zip(*self.truncate_models(std_models,
                                                                                                        truncation_pruning=truncation_pruning,
                                                                                                        truncation_method=truncation_method,
                                                                                                        percent_truncation=percent_truncation,
                                                                                                        residue_scores=zipped_scores,
                                                                                                        work_dir=score_truncate_dir)):
                pre_ensemble = truncated_model[0]
                pre_ensemble_data = copy.copy(truncated_model_data)
                pre_ensemble_data['truncation_score_key'] = score_key.lower()
                
                for ensemble, ensemble_data in zip(*self.edit_side_chains(pre_ensemble,
                                                                          pre_ensemble_data,
                                                                          side_chain_treatments,
                                                                          self.ensembles_directory,
                                                                          single_structure=True)):
                    self.ensembles.append(ensemble)
                    self.ensembles_data.append(ensemble_data)
    
        return self.ensembles
    
    def _generate_residue_scorelist(self, residue_key, score_key, scores):
        """Generate a zipped list of residue indexes and corresponding scores
        
        :residue_key: residue column header keyword
        :score_key: score column header keyword
        :scores: list of dictionaries for each residue
        
        :returns: zipped list of residue index plus score
        """
        
        assert scores[0].has_key(residue_key), "Cannot find residue key in scoresfile"
        assert scores[0].has_key(score_key), "Cannot find score key in scoresfile"
        return [(i[residue_key], i[score_key]) for i in scores]
    
    def _read_scorefile(self, scorefile):
        """
        :scorefile: CSV score file INCLUDING header line
        
        :returns: list of per residue dictionaries containing column data
        """
        scores = []
        with open(scorefile, 'r') as csvfile:
            reader = csv.reader(csvfile)
            header = reader.next()
            # Require this to make sure we have scores matching headers
            for row in reader:
                res = {}
                for h, e in zip(header, row):
                    res[h.lower()] = float(e)
                scores.append(res)
        return scores
