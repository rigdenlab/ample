
from ample.parsers import _contactfile

class EPCMapContactParser(_contactfile.ContactfileParser):
    """ Parser class for EPC-Map contact prediction files """

    _RES1 = 0
    _RES2 = 1
    _RAW_SCORE = 4

    _METHOD = "epcmap"

    def __init__(self):
        _contactfile.ContactfileParser.__init__(self)

    def read(self, contactfile):
        with open(contactfile, 'r') as fh:
            for line in iter(fh.readline, ''): 
                line = line.strip().split()
                
                # Define the contact in a dictionary - use parent method
                contact = self.defineContact(line,
                                             res1_idx=self._RES1,
                                             res2_idx=self._RES2,
                                             raw_score_idx=self._RAW_SCORE,
                                             method=self._METHOD,
                                             file=contactfile)
                self.contacts.append(contact)
        return
