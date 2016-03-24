
import numpy
from ample.parsers import _contactfile_parser

class CCMpredContactParser(_contactfile_parser.ContactfileParser):
    """ Class to parse a CCMpred contact matrix """
    
    _METHOD = "CCMpred"
    
    def __init__(self):
        _contactfile_parser.ContactfileParser.__init__(self)
        
    def read(self, contactmatrix):
        mat = numpy.loadtxt(contactmatrix)

        raw_contacts = self.get_contact_pairs(mat)
        
        for i, j, raw in zip(raw_contacts[0], raw_contacts[1], mat[raw_contacts]):        
            contact = self.contact.copy()
            contact['res1_index'] = i+1             # Matrix starts with 0,
            contact['res2_index'] = j+1             # residue 1 in sequence
            contact['raw_score'] = raw
            contact['method'] = self._METHOD
            contact['file'] = contactmatrix
            self.contacts.append(contact)
        
        self.contacts = self.filter_duplicates(self.contacts)
        
        return
        
    def get_contact_pairs(self, mat):
        """Get the top-scoring contacts"""

        contacts = mat.argsort(axis=None)[::-1]
        contacts = (contacts % mat.shape[0]).astype(numpy.uint16), \
                    numpy.floor(contacts / mat.shape[0]).astype(numpy.uint16)
        return contacts
    
    def filter_duplicates(self, contacts_duplicated):
        
        contacts_filtered = [contact for contact in contacts_duplicated \
                             if contact['res1_index'] < contact['res2_index']]

        return contacts_filtered        


