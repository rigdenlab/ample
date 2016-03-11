"""Test functions for parsers.ccmpred"""
import unittest
from ample.parsers import ccmpred

class Test(unittest.TestCase):
    def test_filter(self):
        c = ccmpred.CCMpredContactParser()
        contacts_duplicated = [{'res1_index': 11, 'res2_index': 100, 'raw_score': 0.5},
                               {'res1_index': 100, 'res2_index': 11, 'raw_score': 0.5},
                               
                               {'res1_index': 20, 'res2_index': 150, 'raw_score': 0.3},
                               {'res1_index': 150, 'res2_index': 20, 'raw_score': 0.3},
                               
                               {'res1_index': 6, 'res2_index': 70, 'raw_score': 0.2},
                               {'res1_index': 70, 'res2_index': 6, 'raw_score': 0.9},
                               
                               {'res1_index': 1, 'res2_index': 8, 'raw_score': 0.2},
                               {'res1_index': 2, 'res2_index': 9, 'raw_score': 0.2},
                               
                               {'res1_index': 50, 'res2_index': 80, 'raw_score': 0.5},
                               
                               {'res1_index': 1, 'res2_index': 1, 'raw_score': -0.000000}]
        
        ref_contacts = [{'res1_index': 11, 'res2_index': 100, 'raw_score': 0.5},
                        {'res1_index': 20, 'res2_index': 150, 'raw_score': 0.3},
                        {'res1_index': 6, 'res2_index': 70, 'raw_score': 0.2},
                        {'res1_index': 1, 'res2_index': 8, 'raw_score': 0.2},
                        {'res1_index': 2, 'res2_index': 9, 'raw_score': 0.2},
                        {'res1_index': 50, 'res2_index': 80, 'raw_score': 0.5}]
        
        
        contacts_filtered = c.filterDuplicates(contacts_duplicated)
        
        self.assertItemsEqual(ref_contacts, contacts_filtered)
    
if __name__ == "__main__":
    unittest.main()