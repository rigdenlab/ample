"""Test functions for ensembler.single_model"""

import unittest
from ample.ensembler import single_model


class Test(unittest.TestCase):
    
    def test_generateResidueScorelist(self):
        ensembler = single_model.Ensembler()
        scores = [{'residue': (i), 'concoord': (i*0.5*3.452), 'rosetta': (i*3.452),
                   'unknown': (i*1.11*3.452)} for i in xrange(1, 11)]
        zipped_concoord = ensembler._generate_residue_scorelist('residue', 'concoord', scores)
        ref_concoord = [(1, 1.726), (2, 3.452), (3, 5.178), (4, 6.904),
                        (5, 8.629999999999999), (6, 10.356), (7, 12.082),
                        (8, 13.808), (9, 15.533999999999999), 
                        (10, 17.259999999999998)]
        self.assertEqual(ref_concoord, zipped_concoord)
        zipped_rosetta = ensembler._generate_residue_scorelist('residue', 'rosetta', scores)
        ref_rosetta = [(1, 3.452), (2, 6.904), (3, 10.356), (4, 13.808),
                       (5, 17.259999999999998), (6, 20.712), (7, 24.164),
                       (8, 27.616), (9, 31.067999999999998), 
                       (10, 34.519999999999996)]
        self.assertEqual(ref_rosetta, zipped_rosetta)
 
if __name__ == "__main__":
    unittest.main()
