"""Test functions for parsers._contactfile_parser"""
import unittest
from ample.parsers import _contactfile_parser

class Test(unittest.TestCase):
    def setUp(self):
        self.cp = _contactfile_parser.ContactfileParser()

    def test_sorting(self):
        contacts = [{'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                    {'res1': 1, 'res2': 100, 'raw_score': -100},
                    {'res1': 6, 'res2': 10, 'raw_score': 1}]

        ref_contacts1 = [{'res1': 1, 'res2': 100, 'raw_score': -100},
                         {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                         {'res1': 6, 'res2': 10, 'raw_score': 1}]

        ref_contacts2 = [{'res1': 6, 'res2': 10, 'raw_score': 1},
                         {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                         {'res1': 1, 'res2': 100, 'raw_score': -100}]

        ref_contacts3 = [{'res1': 6, 'res2': 10, 'raw_score': 1},
                         {'res1': 5, 'res2': 37, 'raw_score': -4.767827},
                         {'res1': 1, 'res2': 100, 'raw_score': -100}]

        self.cp.setContacts(contacts)

        self.cp.sortContacts('res1', False)
        sorted_contacts1 = self.cp.getContacts()
        self.assertEqual(ref_contacts1, sorted_contacts1)

        self.cp.sortContacts('res2', False)
        sorted_contacts2 = self.cp.getContacts()
        self.assertEqual(ref_contacts2, sorted_contacts2)

        self.cp.sortContacts('raw_score', True)
        sorted_contacts3 = self.cp.getContacts()
        self.assertEqual(ref_contacts3, sorted_contacts3)
        
    def test_scalarScoring(self):
        contacts = [{'res1': 5, 'res2': 37, 'raw_score': 0.4},
                    {'res1': 1, 'res2': 100, 'raw_score': 0.6},
                    {'res1': 6, 'res2': 10, 'raw_score': 0.5}]
        
        ref_contacts = [{'res1': 5, 'res2': 37, 'raw_score': 0.4, 'scalar_score': 0.8},
                        {'res1': 1, 'res2': 100, 'raw_score': 0.6 ,'scalar_score': 1.2},
                        {'res1': 6, 'res2': 10, 'raw_score': 0.5, 'scalar_score': 1.0}]
        
        self.cp.setContacts(contacts)
        self.cp.calculateScalarScores()
        scaled_contacts = self.cp.getContacts()
        self.assertEqual(ref_contacts, scaled_contacts)
    
if __name__ == "__main__":
    unittest.main()
