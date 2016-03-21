
from ample.parsers import _contactfile_parser

class BCLContactsContactParser(_contactfile_parser.ContactfileParser):
    """ Parser class for BCL::Contact contact prediction files """
    
    _RES1       = 0
    _RES2       = 2
    _HH_IDX     = 4     # Helix-Helix
    _HS_IDX     = 5     # Helix-Strand
    _SH_IDX     = 6     # Strand-Helix
    _SS_IDX     = 7     # Strand-Strand
    _ShSh_IDX   = 8     # Sheet-Sheet  
    _RAW_SCORE  = 9
    
    _METHOD     = "blccontact"

    def __init__(self):
        _contactfile_parser.ContactfileParser.__init__(self)

    def read(self, contactfile):
        with open(contactfile, 'r') as fh:
            for line in iter(fh.readline, ''):
                if line.startswith("CHAINS"): 
                    continue
                line = line.strip().split()
                
                # Define the contact in a dictionary - use parent method
                contact = self.define_contact(line, 
                                              res1_idx=self._RES1, 
                                              res2_idx=self._RES2, 
                                              raw_score_idx=self._RAW_SCORE,
                                              method=self._METHOD, 
                                              file=contactfile)
                contact['helix_helix'] = line[self._HH_IDX]
                contact['helix_strand'] = line[self._HS_IDX]
                contact['strand_helix'] = line[self._SH_IDX]
                contact['strand_strand'] = line[self._SS_IDX]
                contact['sheet_sheet'] = line[self._ShSh_IDX]
                
                self.contacts.append(contact)
