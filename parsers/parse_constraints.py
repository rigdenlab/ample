#!/usr/bin/env ccp4-python

import parse_contactfile

class ConstraintfileParser(parse_contactfile.ContactfileParser):
    
    def __init__(self):
        parse_contactfile.ContactfileParser.__init__(self)
        
    def read(self, constraintfile):
        """Read a constraints file to convert it to an array of
        contact dictionaries"""

        with open(constraintfile, 'r') as fh:
            for line in iter(fh.readline, ''):

                if line.startswith("AtomPair"):
                    contact = self._atompair(line)
                    contact['method'] = "ample"
                    contact['file'] = constraintfile
                    self.contacts.append(contact)
                
        if not self.contacts:
            msg = "Could not convert constraints to contacts. Unrecognised format."
            raise RuntimeError(msg)
        
        return
    
    def _atompair(self, line):
        """AtomPair specific line extractor"""

        atm1, res1_index, atm2, res2_index = line.strip().split()[1:5]
        
        contact = self.contact.copy()
        contact['atom1'] = atm1
        contact['atom2'] = atm2
        contact['res1_index'] = int(res1_index)
        contact['res2_index'] = int(res2_index)

        return contact
