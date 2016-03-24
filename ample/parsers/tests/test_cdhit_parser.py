"""Test functions for parsers.cdhit_parser"""

import unittest
from ample.parsers import cdhit_parser

class Test(unittest.TestCase):
    def test_average(self):
        data = {0: {'count': 1, 'averageId': 100.00, 'averageLength': 99, 'centroid': 'test1', 'neffWeight': 1.00},
                1: {'count': 3, 'averageId': 332.16, 'averageLength': 290, 'centroid': 'test2', 'neffWeight': 1.00},
                2: {'count': 2, 'averageId': 100.00, 'averageLength': 199, 'centroid': 'test3', 'neffWeight': 1.00},
                3: {'count': 2, 'averageId': 100.00, 'averageLength': 200, 'centroid': 'test4', 'neffWeight': 1.00}}
        
        ref_data = {0: {'count': 1, 'averageId': 100.00, 'averageLength': 99, 'centroid': 'test1', 'neffWeight': 1.00},
                    1: {'count': 3, 'averageId': 110.72, 'averageLength': 96, 'centroid': 'test2', 'neffWeight': 0.3333333333},
                    2: {'count': 2, 'averageId': 50.00, 'averageLength': 99, 'centroid': 'test3', 'neffWeight': 0.5},
                    3: {'count': 2, 'averageId': 50.00, 'averageLength': 100, 'centroid': 'test4', 'neffweight': 0.5}}
        
        cp = cdhit_parser.CDhitLogParser()
        out_data = cp.averageClusters(data)

        self.assertItemsEqual(ref_data, out_data)

if __name__ == "__main__":
    unittest.main()
