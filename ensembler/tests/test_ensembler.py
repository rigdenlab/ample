"""Test functions for ensembler._ensembler"""

import unittest

from ample.ensembler import _ensembler

class Test(unittest.TestCase):
    
    def test_convert_residue_scores(self):
        residue_scores = [(i, 0.1*i) for i in xrange(1, 11)]
        ensembler = _ensembler.Ensembler()
        scores = ensembler._convert_residue_scores(residue_scores)
        score_idxs = [i.idx for i in scores]
        ref_score_idxs = [i for i in xrange(10)] # Minus one compared to org data
        self.assertEqual(ref_score_idxs, score_idxs)
        score_resSeq = [i.resSeq for i in scores]
        ref_score_resSeq = [i for i in xrange(1, 11)] # Same i as org data
        self.assertEqual(ref_score_resSeq, score_resSeq)
        score_variances = [i.variance for i in scores]
        ref_score_variances = [(0.1*i) for i in xrange(1, 11)] # Same i as org data
        self.assertEqual(ref_score_variances, score_variances)