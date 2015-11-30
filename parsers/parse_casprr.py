#!/usr/bin/env ccp4-python

import os

import parse_contactfile


class CaspContactParser(parse_contactfile.ContactfileParser):
    """ Parser class for CASP RR contact prediction files """

    def __init__(self):
        parse_contactfile.ContactfileParser.__init__(self)
        self.method = "CASP RR"

    def checkFormat(self, contactfile):
        """Check for correctness of CASP RR contact file
        :returns: True or False
        """
        with open(contactfile, 'r') as fh: line = fh.readline().strip().split()
        return True if line[0]=="PFRMAT" and line[1]=="RR" else False
        
    def read(self, contactfile):
        self.infile = os.path.abspath(contactfile)
        with open(self.infile, 'r') as fh: self._read(fh)
        return

    def _read(self, fh):
        for line in iter(fh.readline, ''):
            line = line.strip().split()
            # Set method if provided
            if line[0].upper() == "METHOD": self.method = line[1]

            # Ignore headers
            elif line[0].isalpha(): continue

            # Read the contacts
            elif line[0].isdigit():
                
                # Define the contact in a dictionary - use parent method
                contact = self.defineContact(line,
                                             res1_idx=0,
                                             res2_idx=1,
                                             raw_score_idx=4,
                                             method=self.method,
                                             file=self.infile)
                
                contact['lb'] = float(line[2])
                contact['ub'] = float(line[3])
                
                self.contacts.append(contact)
        return
    
    def write(self, outfile):
        assert self.contacts, "No contacts provided"
        
        # Sort the contacts if not done already
        if not self.isSorted: self.sortContacts('res1_index', descending=False)
        
        # Set the name of the method used to get these contacts
        self.method = self.contacts[0]['method']

        # Populate the final output string and write it to a file
        out_str = self._write(self.contacts)
        with open(outfile, 'w') as fh: fh.write(out_str)

        return

    def _write(self, contacts):
        format_str = "PFRMAT RR"
        method_str = "METHOD %s" % self.method.upper()
        model_str  = "MODEL 1"
        # Format the contacts to writable string
        contact_str= "\n".join(self._format_contacts(contacts))
        end_str    = "END"
        out_str = "\n".join([format_str, method_str, model_str, contact_str, end_str])
        return out_str

    def _format_contacts(self, contacts):
        formatted_contacts = [("%4d %4d %2d %2d %f" % (contact['res1_index'],
                                                       contact['res2_index'],
                                                       contact['lb'], 
                                                       contact['ub'],
                                                       contact['raw_score'])) \
                                    for contact in contacts]
        return formatted_contacts
