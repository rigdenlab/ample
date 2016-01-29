#!/usr/bin/env ccp4-python

import parse_contactfile
      

class RestraintfileParser(parse_contactfile.ContactfileParser):
    """Parser class for Rosetta restraints files"""
    
    def __init__(self):
        parse_contactfile.ContactfileParser.__init__(self)
        
    def read(self, restraintfile):
        """Read a restraints file to convert it to an array of
        contact dictionaries"""

        with open(restraintfile, 'r') as fh:
            for line in iter(fh.readline, ''):
                contact = self.contact.copy()       # Universal contact template

                # Read in the contact data 
                if not line.strip() or line.startswith("#"):
                    continue
                elif line.startswith("AtomPair"): 
                    self._atompair(contact, line)   # AtomPair specific function
                else: 
                    msg = "Unrecognised format: {0}".format(line.strip())
                    raise ValueError(msg)
                
                # Fulfill with remaining metadata information
                contact['method'] = "ample"
                contact['file'] = restraintfile
                self.contacts.append(contact)       # Store contact
                
        if not self.contacts:
            msg = "No converted restraints"
            raise ValueError(msg)
        
        return
    
    def _atompair(self, contact, line):
        """AtomPair specific line extractor"""

        atm1, res1_index, atm2, res2_index = line.strip().split()[1:5]
        
        contact['atom1'] = atm1
        contact['atom2'] = atm2
        contact['res1_index'] = int(res1_index)
        contact['res2_index'] = int(res2_index)

        return
