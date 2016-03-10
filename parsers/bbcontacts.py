
from ample.parsers import _contactfile

class BBcontactsContactParser(_contactfile.ContactfileParser):
    """ Parser class for bbcontacts contact prediction files """

    _ID         = 0
    _DIV_FACTOR = 1
    _STRAND_ORI = 2
    _RAW_SCORE = 3
    _STRAND_IDX = 4
    _STRAND_POS = 5
    _RES1       = 7     # Annoyingly are the other way around in file
    _RES2       = 6     # So we swap them here for aesthetic purposes
    
    _METHOD     = "bbcontacts"
    
    def __init__(self):
        _contactfile.ContactfileParser.__init__(self)

    def read(self, contactfile, removeFirstLast=True):
        with open(contactfile, 'r') as fh:
            for line in iter(fh.readline, ''):
                if line.startswith("#") or any(i.upper()=="NA" for i in line.split()): 
                    continue
                line = line.strip().split()
               
                ## Remove two-residue strands. High FP rate
                if line[self._STRAND_POS].upper() != "FIRST" \
                    and removeFirstLast \
                    and self.isFirstLast(line):
                        self.contacts.pop()
                        continue
                
                # Define the contact in a dictionary - use parent method
                contact = self.defineContact(line, 
                                             res1_idx=self._RES1, 
                                             res2_idx=self._RES2, 
                                             raw_score_idx=self._RAW_SCORE,
                                             method=self._METHOD, 
                                             file=contactfile)
                
                # Additional bbcontacts specific attributes
                contact['diversity_factor'] = float(line[self._DIV_FACTOR])
                contact['strand_orientation'] = line[self._STRAND_ORI]
                contact['strand_index'] = int(line[self._STRAND_IDX])
                contact['internal_strand_position'] = line[self._STRAND_POS]
            
                self.contacts.append(contact)
        return

    def isFirstLast(self, line):
        return True \
            if self.contacts[-1]['internal_strand_position'].upper() == "FIRST" \
            and line[self._STRAND_POS].upper() == "LAST" \
            else False
##End BBcontactsContactParser

