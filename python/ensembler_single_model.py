#!/usr/bin/env ccp4-python

"""
16.02.2016

@author: hlfsimko
"""

# System
import copy
import csv
import logging
import sys
import os
import unittest

# Custom
import ample_util
import ensembler
import pdb_edit

# Inherit some other variables defined in the ensembler module
SIDE_CHAIN_TREATMENTS = ensembler.SIDE_CHAIN_TREATMENTS

_logger = logging.getLogger(__name__)


class Ensembler(ensembler.Ensembler):
    """Ensemble creator using on a single input structure and a corresponding
       score file with per residue scores for truncation
    """
        
    def __init__(self):
        # Inherit all functions from Parent Ensembler
        ensembler.Ensembler.__init__(self)
        
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
            raise RuntimeError("More than 1 structure provided")
        
        if len(truncation_scorefile_header) < 2:
            raise RuntimeError("At least two header options for scorefile are required")

        # standardise the structure
        std_models_dir = os.path.join(work_dir, "std_models")
        os.mkdir(std_models_dir)
        std_models = []
        for m in models:
            std_model = ample_util.filename_append(m, 'std', std_models_dir)
            pdb_edit.standardise(pdbin=m, pdbout=std_model, del_hetatm=True)
            std_models.append(std_model)
        _logger.info('Standardised input model: {0}'.format(std_models[0]))
              
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)

        truncate_dir = os.path.join(self.work_dir,"single_truncate")
        if not os.path.isdir(truncate_dir): os.mkdir(truncate_dir)
        
        # Read all the scores into a per residue dictionary
        residue_scores = self._read_scorefile(truncation_scorefile)
        residue_key = truncation_scorefile_header.pop(0).lower()
        
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


class Test(unittest.TestCase):
    
    def test_generateResidueScorelist(self):
        """Test the extraction of one score per residue"""
        ########################################################################
        # Input data
        ########################################################################    
        scores = [{'residue': (i), 
                   'concoord': (i*0.5*3.452), 
                   'rosetta': (i*3.452),
                   'unknown': (i*1.11*3.452)} \
                    for i in xrange(1, 11)]
        
        ########################################################################
        # Reference data
        ########################################################################
        ref_concoord = [(1, 1.726), (2, 3.452), (3, 5.178), (4, 6.904),
                        (5, 8.629999999999999), (6, 10.356), (7, 12.082),
                        (8, 13.808), (9, 15.533999999999999), 
                        (10, 17.259999999999998)]
        
        ref_rosetta = [(1, 3.452), (2, 6.904), (3, 10.356), (4, 13.808),
                       (5, 17.259999999999998), (6, 20.712), (7, 24.164),
                       (8, 27.616), (9, 31.067999999999998), 
                       (10, 34.519999999999996)]
        
        ########################################################################
        # Function testing
        ######################################################################## 
        ensembler = Ensembler()
        
        zipped_concoord = ensembler._generate_residue_scorelist('residue', 'concoord', scores)
        zipped_rosetta = ensembler._generate_residue_scorelist('residue', 'rosetta', scores)
        
        self.assertEqual(ref_concoord, zipped_concoord)
        self.assertEqual(ref_rosetta, zipped_rosetta)
        
        return
        

if __name__ == "__main__":
    import argparse
    
    options = argparse.ArgumentParser()
    options.add_argument('-m', type=str, nargs=1, dest='model', required=True,
                         help='path to model')
    options.add_argument('-n', type=int, dest='nproc',
                         help='number of processors')
    options.add_argument('-p', type=float, 
                         dest='percent_truncation', default=5, 
                         help='truncation percentage')
    options.add_argument('-s', type=str, required=True,
                         dest='truncation_scorefile',
                         help='truncation scorefile')
    options.add_argument('-sh', nargs='+', required=True,
                         dest='truncation_scorefile_header',
                         help='truncation scorefile header')
    options.add_argument('-w', type=str, dest='work_dir',
                         help='working directory')
    args = vars(options.parse_args())
    
    args['truncation_method'] = 'scores'
    
    models = [ os.path.abspath(i) for i in args['model'] ]
    args['truncation_scorefile'] = os.path.abspath(args['truncation_scorefile'])
    
    if not args['work_dir']: args['work_dir'] = os.getcwd()
    
    ensembles_directory = os.path.join(args['work_dir'], 'ensembles')
    work_dir = os.path.join(args['work_dir'], "ensemble_workdir")
    
    if not os.path.isdir(ensembles_directory): os.mkdir(ensembles_directory)
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    os.chdir(work_dir)

    ensembler = Ensembler()
    ensembles = ensembler.generate_ensembles(models,
                                             percent_truncation=args['percent_truncation'],
                                             truncation_method=args['truncation_method'],
                                             truncation_scorefile=args['truncation_scorefile'],
                                             truncation_scorefile_header=args['truncation_scorefile_header'],
                                             nproc=args['nproc'],
                                             ensembles_directory=ensembles_directory,
                                             work_dir=work_dir)
      
    args['ensemble_ok'] = os.path.join(args['work_dir'], "ensemble.ok")
    with open(args['ensemble_ok'],'w') as f: f.write('ok\n')
    
    
    
