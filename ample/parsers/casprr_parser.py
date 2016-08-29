
import os

from ample.parsers import _contactfile_parser

__author__ = "Felix Simkovic"
__date__ = "16.01.2016"
__version__ = "1.0"

class CaspContactParser(_contactfile_parser.ContactfileParser):
    """Parser class for CASP RR contact prediction files
    
    This class servers as a parser for contact prediction files in 
    CASP RR format. It allows you to read and write contact files
    as well as check their format.
    """

    def __init__(self):
        _contactfile_parser.ContactfileParser.__init__(self)
        self.method = "CASP RR"

    def check_format(self, contactfile):
        """This function checks the format of a contact prediction file
    
        This function can be used to check if a contact prediction file 
        is in CASP RR format based on the first line. This must correspond 
        to ``PFRMAT RR``.

        Args:
            contactfile (str): The path to a contact file to read

        Returns:
            bool: True if file is in CASP RR format, else False
        """
        with open(contactfile, 'r') as fh: line = fh.readline().strip().split()
        if line[0]=="PFRMAT" and line[1]=="RR":
            return True
        return False
        
    def read(self, contactfile):
        """Read a contactfile into a list of contacts
        
        Args:
            contactfile (str): The path to a contact file to read
        """
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
                contact = self.define_contact(line,
                                             res1_idx=0,
                                             res2_idx=1,
                                             raw_score_idx=4,
                                             method=self.method,
                                             file=self.infile)
                
                contact['lb'] = float(line[2])
                contact['ub'] = float(line[3])
                
                self.contacts.append(contact)
        return
    
    def write(self, outfile, score="raw_score"):
        """Write a set of contacts to a file
        
        Write a list of contacts to a file, sorted by a the score provided.

        Args:
            outfile (str): The path to the CASP RR output file
            score (str): The value by which the contacts will be sorted
        """

        assert self.contacts, "No contacts provided"
        # Sort the contacts if not done already
        if not self.isSorted: self.sort_contacts('res1_index', descending=False)
        # Set the name of the method used to get these contacts
        self.method = self.contacts[0]['method']
        # Populate the final output string and write it to a file
        out_str = self._write(self.contacts, score)
        with open(outfile, 'w') as fh: fh.write(out_str)

        return

    def _write(self, contacts, score):
        format_str = "PFRMAT RR"
        method_str = "METHOD %s" % self.method.upper()
        model_str  = "MODEL 1"
        # Format the contacts to writable string
        contact_str= "\n".join(self._format_contacts(contacts, score))
        end_str    = "END"
        out_str = "\n".join([format_str, method_str, model_str, contact_str, end_str])
        return out_str

    def _format_contacts(self, contacts, score):
        formatted_contacts = [("%4d %4d %2d %2d %f" % (contact['res1_index'],
                                                       contact['res2_index'],
                                                       contact['lb'], 
                                                       contact['ub'],
                                                       contact[score])) \
                                    for contact in contacts]
        return formatted_contacts
