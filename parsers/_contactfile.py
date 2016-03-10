
import operator

class ContactfileParser(object):
    def __init__(self):
        self.contacts = []
        self.infile = None
        self.outfile = None
        
        self.isSorted = False
        
        self.contact = self._contactTemplate()

    def assignAminoAcids(self, sequence):
        assert self.contacts, "No contacts defined"
        
        for contact in self.contacts:
            # Assign the amino acids to each contact
            contact['res1'] = sequence[ contact['res1_index']-1 ]
            contact['res2'] = sequence[ contact['res2_index']-1 ]
            
            # Whilst doing that we can also do the atoms, i.e. Ca or Cb
            contact['atom1'] = "CA" if contact['res1']=="G" else "CB"
            contact['atom2'] = "CA" if contact['res2']=="G" else "CB"
            
        return

    def calculateScalarScores(self):
        """Calculate the scale score for each contact
    
        Each scale score is defined by (raw_score/average(raw_scores))
        """
        raw_scores = [i['raw_score'] for i in self.contacts]
        avgRawScore = sum(raw_scores) / len(raw_scores)
        
        for contact in self.contacts:
            raw_score = float(contact['raw_score'])
            contact['scalar_score'] = raw_score / avgRawScore
        
        return
    
    def defineContact(self, line, **kwargs):
        contact = self.contact.copy()
        
        try:
            contact['res1_index'] = int(line[kwargs['res1_idx']])
            contact['res2_index'] = int(line[kwargs['res2_idx']])
            contact['raw_score'] = float(line[kwargs['raw_score_idx']])
            contact['method'] = kwargs['method']
            contact['file'] = kwargs['file']
        except (ValueError, IndexError):
            msg = "Cannot process contact, please check file format"
            raise RuntimeError(msg)
        
        return contact

    def getContacts(self):
        return self.contacts

    def setContacts(self, contact_list):
        self.contacts = contact_list

    def sortContacts(self, key, descending=False):
        assert self.contacts, "No contacts provided"
        assert key in self.contacts[0], "Key not defined"

        self.contacts = sorted(self.contacts,
                               key=operator.itemgetter(key),
                               reverse=descending)
        self.isSorted = True
        return
    
    def _contactTemplate(self):
        """ create a contact template """
        
        d = {"atom1": None,
             "atom2": None,
             "raw_score": 0.0,
             "scalar_score": 0.0,
             "diversity_factor": 0.0,
             "file": None,
             "internal_strand_position": None,
             "lb": 0.,
             "method": None,
             "res1": None,
             "res2": None,
             "res1_index": 0,
             "res2_index": 0,
             "strand_index": 0,
             "strand_orientation": None,
             "true_positive": True,
             "ub": 8.,
             "weight": 1
            }
        
        return d


