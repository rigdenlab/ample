#!/usr/bin/env ccp4-python

import operator
import unittest

class ContactfileParser(object):
    def __init__(self):
        self.contacts = []
        self.infile = None
        self.outfile = None
        
        self.isSorted = False
        
        self.contact = self._contactTemplate()

    def assignAminoAcids(self, sequence):
        assert self.contacts, "No contacts defined"
        
        for contact in self.contacts:
            # Assign the amino acids to each contact
            contact['res1'] = sequence[ contact['res1_index']-1 ]
            contact['res2'] = sequence[ contact['res2_index']-1 ]
            
            # Whilst doing that we can also do the atoms, i.e. Ca or Cb
            contact['atom1'] = "CA" if contact['res1']=="G" else "CB"
            contact['atom2'] = "CA" if contact['res2']=="G" else "CB"
            
        return

    def calculateScalarScores(self):
        """Calculate the scale score for each contact
    
        Each scale score is defined by (raw_score/average(raw_scores))
        """
        raw_scores = [i['raw_score'] for i in self.contacts]
        avgRawScore = sum(raw_scores) / len(raw_scores)
        
        for contact in self.contacts:
            raw_score = float(contact['raw_score'])
            contact['scalar_score'] = raw_score / avgRawScore
        
        return
    
    def defineContact(self, line, **kwargs):
        contact = self.contact.copy()
        
        try:
            contact['res1_index'] = int(line[kwargs['res1_idx']])
            contact['res2_index'] = int(line[kwargs['res2_idx']])
            contact['raw_score'] = float(line[kwargs['raw_score_idx']])
            contact['method'] = kwargs['method']
            contact['file'] = kwargs['file']
        except (ValueError, IndexError):
            msg = "Cannot process contact, please check file format"
            raise RuntimeError(msg)
        
        return contact

    def getContacts(self):
        return self.contacts

    def setContacts(self, contact_list):
        self.contacts = contact_list

    def sortContacts(self, key, descending=False):
        assert self.contacts, "No contacts provided"
        assert key in self.contacts[0], "Key not defined"

        self.contacts = sorted(self.contacts,
                               key=operator.itemgetter(key),
                               reverse=descending)
        self.isSorted = True
        return
    
    def _contactTemplate(self):
        """ create a contact template """
        
        d = {"atom1": None,
             "atom2": None,
             "raw_score": 0.0,
             "scalar_score": 0.0,
             "diversity_factor": 0.0,
             "file": None,
             "internal_strand_position": None,
             "lb": 0.,
             "method": None,
             "res1": None,
             "res2": None,
             "res1_index": 0,
             "res2_index": 0,
             "strand_index": 0,
             "strand_orientation": None,
             "true_positive": True,
             "ub": 8.,
             "weight": 1
            }
        
        return d


class Test(unittest.TestCase):
    def setUp(self):
        self.cp = ContactfileParser()

    def testSorting(self):
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
        
    def testScalarScoring(self):
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

