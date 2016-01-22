#!/usr/bin/env ccp4-python

import numpy
import parse_contactfile
import unittest


class CCMpredContactParser(parse_contactfile.ContactfileParser):
    """ Class to parse a CCMpred contact matrix """
    
    _METHOD = "CCMpred"
    
    def __init__(self):
        parse_contactfile.ContactfileParser.__init__(self)
        
    def read(self, contactmatrix):
        mat = numpy.loadtxt(contactmatrix)

        raw_contacts = self.getContactPairs(mat)
        
        for i,j, raw in zip(raw_contacts[0], raw_contacts[1], mat[raw_contacts]):        
            contact = self.contact.copy()
            contact['res1_index'] = i+1             # Matrix starts with 0,
            contact['res2_index'] = j+1             # residue 1 in sequence
            contact['raw_score'] = raw
            contact['method'] = self._METHOD
            contact['file'] = contactmatrix
            self.contacts.append(contact)
        
        self.contacts = self.filterDuplicates(self.contacts)
        
        return
        
    def getContactPairs(self, mat):
        """Get the top-scoring contacts"""

        contacts = mat.argsort(axis=None)[::-1]
        contacts = (contacts % mat.shape[0]).astype(numpy.uint16), \
                    numpy.floor(contacts / mat.shape[0]).astype(numpy.uint16)
        return contacts
    
    def filterDuplicates(self, contacts_duplicated):
        
        contacts_filtered = [contact for contact in contacts_duplicated \
                             if contact['res1_index'] < contact['res2_index']]

        return contacts_filtered        


class Test(unittest.TestCase):
    def testFilter(self):
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
        
        c = CCMpredContactParser()
        contacts_filtered = c.filterDuplicates(contacts_duplicated)
        
        self.assertItemsEqual(ref_contacts, contacts_filtered)
    
