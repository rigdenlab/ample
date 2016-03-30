"""Test functions for parsers._contactfile_parser"""

import unittest
from ample.parsers import _contactfile_parser

class Test(unittest.TestCase):

    def test_sorting(self):
        contacts = [{'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                    {'res1': 1, 'res2': 100, 'raw_score': -100},
                    {'res1': 6, 'res2': 10, 'raw_score': 1}]
        cp = _contactfile_parser.ContactfileParser()
        
        cp.set_contacts(contacts)
        cp.sort_contacts('res1', False)
        sorted_contacts = cp.get_contacts()
        ref_contacts = [{'res1': 1, 'res2': 100, 'raw_score': -100},
                        {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                        {'res1': 6, 'res2': 10, 'raw_score': 1}]
        self.assertEqual(ref_contacts, sorted_contacts)
        
        cp.set_contacts(contacts)
        cp.sort_contacts('res2', False)
        sorted_contacts = cp.get_contacts()
        ref_contacts = [{'res1': 6, 'res2': 10, 'raw_score': 1},
                        {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                        {'res1': 1, 'res2': 100, 'raw_score': -100}]
        self.assertEqual(ref_contacts, sorted_contacts)

        cp.set_contacts(contacts)
        cp.sort_contacts('raw_score', True)
        sorted_contacts = cp.get_contacts()
        ref_contacts = [{'res1': 6, 'res2': 10, 'raw_score': 1},
                        {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                        {'res1': 1, 'res2': 100, 'raw_score': -100}]
        self.assertEqual(ref_contacts, sorted_contacts)
        
    def test_scalar_scoring(self):
        contacts = [{'res1': 5, 'res2': 37, 'raw_score': 0.4},
                    {'res1': 1, 'res2': 100, 'raw_score': 0.6},
                    {'res1': 6, 'res2': 10, 'raw_score': 0.5}]
        cp = _contactfile_parser.ContactfileParser()
        cp.set_contacts(contacts)
        cp.calculate_scalar_scores()
        scaled_contacts = cp.get_contacts()
        ref_contacts = [{'res1': 5, 'res2': 37, 'raw_score': 0.4, 'scalar_score': 0.8},
                        {'res1': 1, 'res2': 100, 'raw_score': 0.6 ,'scalar_score': 1.2},
                        {'res1': 6, 'res2': 10, 'raw_score': 0.5, 'scalar_score': 1.0}]
        self.assertEqual(ref_contacts, scaled_contacts)
    
if __name__ == "__main__":
    unittest.main()
