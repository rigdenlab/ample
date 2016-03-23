
from ample.parsers import _contactfile_parser

class EVfoldContactParser(_contactfile_parser.ContactfileParser):
    """ Parser class for evfold contact prediction files """

    _RES1 = 0
    _RES2 = 2
    _RAW_SCORE = 5

    _METHOD = "evfold"

    def __init__(self):
        _contactfile_parser.ContactfileParser.__init__(self)

    def read(self, contactfile):
        with open(contactfile, 'r') as f:
            for line in iter(f.readline, ''):
                line = line.strip().split()
                if line[0].isdigit():
                    
                    # Define the contact in a dictionary - use parent method
                    contact = self.define_contact(line,
                                                  res1_idx=self._RES1,
                                                  res2_idx=self._RES2,
                                                  raw_score_idx=self._RAW_SCORE,
                                                  method=self._METHOD,
                                                  file=contactfile)
                    
                    self.contacts.append(contact)
        return
