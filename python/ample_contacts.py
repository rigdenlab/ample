#!/usr/bin/env ccp4-python

'''
14.11.2015

@author: hlfsimko
'''

# System
import logging
import numpy
import os
import sys
import unittest

if "CCP4_AMPLE_ROOT" in os.environ.keys() and "CCP4" in os.environ.keys(): 
    root = os.environ["CCP4_AMPLE_ROOT"]
elif "CCP4" in os.environ.keys(): 
    root = os.path.join(os.environ["CCP4"], "share", "ample")
else:
    raise RuntimeError('CCP4 not found')

sys.path.insert(0, os.path.join(root, "parsers"))

# Custom
import ample_sequence
import energy_functions
import parse_casprr
import parse_constraints
import parse_psipred
import pdb_edit

# Check matplotlib is available - problems with Windows
try:
    import ample_plot
    _MATPLOTLIB = True
except ImportError:
    _MATPLOTLIB = False

# Check biopython is availbale - compatibility with older ccp4-python 
try:
    import parse_alignment
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False


def checkOptions(amoptd):
    """Function to check that all contact files are available"""
    
    # Make sure contact file is provided with bbcontacts_file
    if not amoptd['contact_file'] and amoptd['bbcontacts_file']:
        msg = "Must provide -contact_file when using -bbcontacts_file or use as -contact_file"
        raise RuntimeError(msg)
    
    # Check the existence of the contact file and whether it is CASP RR format
    # based on the required header `PFRMAT RR`
    if amoptd['contact_file']:
        if not os.path.exists(str(amoptd['contact_file'])):
            msg = "Cannot find contact file:\n{0}".format(amoptd['contact_file'])
            raise RuntimeError(msg)
           
        if not parse_casprr.CaspContactParser().checkFormat(amoptd['contact_file']):
            msg = "Wrong format in contact file:\n{0}".format(amoptd['contact_file'])
            raise RuntimeError(msg)
    
    # Check the existence of the contact file and whether it is CASP RR format
    # based on the required header `PFRMAT RR`   
    if amoptd['bbcontacts_file']:
        if not os.path.exists(amoptd['bbcontacts_file']):
            msg = "Cannot find contact file:\n{0}".format(amoptd['contact_file'])
            raise RuntimeError(msg)
            
        if not parse_casprr.CaspContactParser().checkFormat(amoptd['bbcontacts_file']):
            msg = "Wrong format in contact file:\n{0}".format(amoptd['bbcontacts_file'])
            raise RuntimeError(msg)
    
    # Make sure user selected energy function is pre-defined
    if amoptd['energy_function']:
        try: 
            energyFunction = getattr(energy_functions, amoptd['energy_function'])
        except AttributeError:
            msg = "Rosetta energy function {0} unavailable".format(amoptd['energy_function'])
            raise RuntimeError(msg)
    
    return


class Contacter(object):
    """ Class to handle contact predictions """
    
    def __init__(self, optd=None):
        self.logger = logging.getLogger()
        
        self.bbcontacts_file = None
        self.constraints_file = None
        self.contact_file = None
        self.contact_map = None
        self.contact_ppv = None
        self.psipred_ss2 = None
        self.structure_pdb = None
        
        self.contacts = None
        self.raw_contacts = None
        self.sequence = None
        
        if optd: self.init(optd)
        
        return
    
    def init(self, optd):
        self.optd = optd

        self.contact_file = optd['contact_file']
        self.contact_map = os.path.join(optd['work_dir'], optd['name'] + ".cm.pdf")
        self.constraints_file = os.path.join(optd['work_dir'], optd['name'] + ".cst") \
            if not optd['constraints_file'] else optd['constraints_file']
        self.sequence = optd['sequence']
        
        # Optional files
        if optd['native_pdb_std']: self.structure_pdb=optd['native_pdb']
        
        
        if optd['psipred_ss2']: self.psipred_ss2=optd['psipred_ss2']
        
        # Extract the raw sequence
        if self.contact_file:
            self.raw_contacts = self._readContacts(self.contact_file, self.sequence)
            # Prepare the contacts depending on all the user/default options
            self.contacts = self._prepare(self.raw_contacts, 
                                          optd['constraints_factor'], 
                                          optd['distance_to_neighbour'])
            
            # Map bbcontacts on top of the previously provided contacts
            if optd['bbcontacts_file']: 
                self.bbcontacts_file = optd['bbcontacts_file']
                self._readAdditionalBBcontacts(optd['bbcontacts_file'])
                
        else:
            self.contacts = self._readConstraints(self.constraints_file)
   
        return
                
    def format(self, constraintfile):
        """ Format contacts to Rosetta constraints """
        
        self.logger.info("Re-formatting contacts to constraints using the {0} function".format(self.optd['energy_function']))
        
        # Format the contacts to constraints
        contact_formatted_lines = self._formatToConstraints(self.contacts, self.optd['energy_function'])
        
        # Write contacts to constraint file
        with open(constraintfile, 'w') as oh: oh.write("\n".join(contact_formatted_lines))
        
        return
       
    def _formatToConstraints(self, contacts, user_function):
        """ Return a list of Rosetta string lines """
    
        try: 
            energyFunction = getattr(energy_functions, user_function)
        except AttributeError:
            msg = "Rosetta energy function `{0}` unavailable".format(user_function)
            raise RuntimeError(msg)
        
        # Format each contact according to the line provided above 
        contact_formatted_lines = [ energyFunction(contact) \
                                        for contact in contacts ]
        
        return contact_formatted_lines
       
    def process_contactfile(self):
        """Wrapper function for
            1) contact formatting
            2) contact map plotting
            3) calculation of contact accuracy (PPV)  
        """
        
        assert self.contacts, "No contacts provided"

        # Format contacts
        self.format(self.constraints_file)
           
        # Contact map plotting. We can parse ss2file and structure file blindly, checks in place
        self.plot(self.contact_map,
                  ss2file=self.psipred_ss2,
                  structurefile=self.structure_pdb)
        
        if self.structure_pdb: 
            self.contact_ppv=self.ppv(self.structure_pdb)
        
        return
    
    def process_constraintsfile(self):
        """Wrapper function for
            1) contact map plotting
            2) calculation of contact accuracy (PPV)
        """
        if not self.contacts: return
        
        # Contact map plotting. We can parse ss2file and structure file blindly, checks in place
        self.plot(self.contact_map,
                  ss2file=self.psipred_ss2,
                  structurefile=self.structure_pdb)
        
        if self.structure_pdb:
            self.contact_ppv=self.ppv(self.structure_pdb)
        
        return
          
    def plot(self, figurefile, ss2file=None, structurefile=None, offset=0):
        """ Plot a contact map """
        
        # make sure we can import matplotlib
        if not _MATPLOTLIB: 
            self.logger.warning("Cannot plot contact map due to missing python dependencies")
            return
        
        # Just to make sure we have a structurefile and we can import Biopython
        availableStructure = True if structurefile and _BIOPYTHON else False
            
        ap = ample_plot.Plotter()
        ap.initialise(figsize=(5,5), dpi=600)
        
        # Get the coordinates from the reference structure and plot in gray
        if availableStructure:
            RCs, RCm = self.structureContacts(structurefile, self.sequence)
            RCs = list([i+offset for i in R] for R in RCs)
            ap.addScatter(RCs[0], RCs[1], marker='.', c='#DDDDDD', s=10, \
                                                      edgecolor="none", linewidths=0.0)
            
        # Get the coordinates for the secondary structure prediction and plot it alogn the diagonal
        if ss2file:
            SCs, Scolors = self.secondaryStructureContacts(ss2file)
            SCs = list([i+offset for i in S] for S in SCs)
            ap.addScatter(SCs[0], SCs[1], marker='.', c=Scolors, s=35, \
                                                      edgecolor="black", linewidths=0.1)
        
        # Get the TP and FP colours
        tp_colors = self._tp_codes(self.contacts, RCm, self.structure_seq, offset=offset) \
                        if availableStructure \
                            else ['#004F9D']

        # Bit cleaner code if we extract X and Y coordinates and then parse them two the plotter
        # Reduce residue index by one to account for counting from 0
        Xs = [ i['res1_index']+offset-1 for i in self.contacts ]
        Ys = [ i['res2_index']+offset-1 for i in self.contacts ]
 
        ap.addScatter(Xs, Ys, marker='.', c=tp_colors, s=10, \
                                          edgecolor="none", linewidths=0.0)
        ap.addScatter(Ys, Xs, marker='.', c=tp_colors, s=10, \
                                          edgecolor="none", linewidths=0.0)
        
        ap.axisLimits(len(self.sequence), offset=offset)
        ap.axisTitles("Residue number", "Residue number")

        ap.saveFig(figurefile)
        
        return
    
    def _tp_codes(self, contacts, RCm, RCm_sequence, offset=0):
        '''Get the color codes for each contact depending on match and weight'''
        
        tp_colors = []
        
        for idx in range(len(contacts)):
            c_x = contacts[idx]['res1_index']-offset-1
            c_y = contacts[idx]['res2_index']-offset-1

            if RCm[c_x, c_y] > 0    and contacts[idx]['weight']==2: tp_colors.append('#2D9D00')
            elif RCm[c_x, c_y] == 0 and contacts[idx]['weight']==2: tp_colors.append('#AB0000')
            elif RCm[c_x, c_y] > 0  and contacts[idx]['weight']==1: tp_colors.append("#38C700")
            elif RCm[c_x, c_y] == 0 and contacts[idx]['weight']==1: tp_colors.append("#D70909")
            else: tp_colors.append('#004F9D')
            
        return tp_colors
    
    def ppv(self, structurefile):
        """Calculate the True Positive Rate of contact prediction"""
         
        # Check that we were able to import Biopython
        if not _BIOPYTHON: return 0.0
        
        assert structurefile, "You need to provide PDB to calculate the PPV"
        
        RCs, RCm = self.structureContacts(structurefile, self.sequence)
        
        ppv = self._ppv_score(self.contacts, RCm, self.structure_seq)
        
        self.logger.info("Accuracy of contact prediction (PPV): {0} %".format(ppv*100))
        
        return ppv
    
    def _ppv_score(self, contacts, RCm, RCm_sequence, offset=0):
        '''Calculate the Positive Predicted Value (PPV) for the predicted contacts'''
        
        GP, FP, TP = 0.0, 0.0, 0.0
        
        for idx in range(len(contacts)):
            c_x = contacts[idx]['res1_index']-offset-1
            c_y = contacts[idx]['res2_index']-offset-1
            if RCm_sequence[c_x] == '-' or RCm_sequence[c_y] == '-': GP += 1.0
            elif RCm[c_x, c_y] > 0: TP += 1.0
            else: FP += 1.0
        
        ppv = TP/(TP+FP)
        
        assert (TP+FP+GP)==len(contacts), "Differing number of contacts used for PPV calculation"
        
        return ppv
    
    def secondaryStructureContacts(self, ss2file):
        pp = parse_psipred.PsipredSs2Parser(ss2file)
        pred = pp.getSecondaryStructure()
        
        coords = [[],[]]
        colors = []
        for i in xrange(len(pred)):
            coords[0].append(i)
            coords[1].append(i)
            if   pred[i]=="H": colors.append('#8B0043')
            elif pred[i]=="E": colors.append('#0080AD')
            elif pred[i]=="C": colors.append('#CCCCCC')
        return coords, colors
    
    def structureContacts(self, structure, alignmentSequence=None):
        """ Extract all contacts from a PDB structure.
            Optional: provide an additional sequence to adjust the structure contacts to.
        """
        
        self.structure_seq = pdb_edit.sequence(structure).values()[0]

        cb_lst = [numpy.array(x) for x in pdb_edit.xyz_cb_coordinates(structure)]
        
        # Adjust the residue list to that of the input sequence
        if alignmentSequence and _BIOPYTHON:
            aligned_seq_list = parse_alignment.AlignmentParser().align_sequences(alignmentSequence, self.structure_seq)
   
            j = 0
            gapped_cb_lst=[]
            for i in xrange(len(aligned_seq_list[1])):
                if aligned_seq_list[1][i] == '-': gapped_cb_lst.append('-')
                elif aligned_seq_list[0][i] == '-': j += 1
                else:
                    gapped_cb_lst.append(cb_lst[j])
                    j += 1
            cb_lst = gapped_cb_lst
            
            self.structure_seq = aligned_seq_list[1]
        
        contacts, contact_map = self._cb_contacts(cb_lst, cb_lst, len(self.structure_seq))

        return contacts, contact_map
    
    def _cb_contacts(self, cb1_lst, cb2_lst, length, cutoff=8):
        '''Get the contacts between the two lists of contacts'''

        self.logger.info("Distance cutoff of participating atoms is: %dA" % cutoff)

        dist_mat = numpy.zeros((length, length), numpy.float)
        dist_mat.fill(float('inf'))
        
        for i, cb1 in enumerate(cb1_lst):
            for j, cb2 in enumerate(cb2_lst):
                
                # Skip infinite against infinite vector. Result of aligned sequence
                #  without any contacts in that region. 
                if cb1 == '-' \
                    or cb2 == '-' \
                    or (numpy.isinf(cb1).all() and numpy.isinf(cb2).all()): 
                        continue
                
                diff_vec = cb1 - cb2
                # Square all to make them positive, get sum and then revert with sqrt
                dist_mat[i,j] = numpy.sqrt(numpy.sum(diff_vec ** 2))

        # Avoid warning messages to the console in this particular area 
        # regarding the `inf` values in the dist_mat matrix. 
        numpy.seterr(invalid='ignore')

        ref_map = dist_mat < cutoff
        ref_contacts = numpy.where(dist_mat < cutoff)
        
        numpy.seterr(invalid='print')    # Reset the warning message board

        return ref_contacts, ref_map        
    
    def _prepare(self, contacts, constraint_factor, distance_to_neighbour):
        # No processing so far, but we will need to do it to match atom to contact
        nrConstraints = int(len(self.sequence) * constraint_factor)
        
        # Filter out any contact pairs that are within `distance residues` of each other
        contacts_filtered = self._filterNeighbours(contacts, distance_to_neighbour)
        
        # Truncate to the number of constraints defined by the user
        contacts_filtered_truncated = self._truncateContactList(contacts_filtered, nrConstraints)
        
        return contacts_filtered_truncated
   
    def _filterNeighbours(self, contacts, distance):
        ''' Filter any contact pairs closer than `distance` residues '''
        contacts_filtered = [contact for contact in contacts \
                                if abs(contact['res2_index'] - contact['res1_index']) >= distance]
        return contacts_filtered  
    
    def _truncateContactList(self, contacts, nrConstraints):
        return contacts[:nrConstraints]
    
    def _readAdditionalBBcontacts(self, secondContactfile):
        """ Method designated solely for reading bbcontacts to overlay/match
            them to the existing set of contacts.
            
            *** BE CAREFUL WHEN USING WITH ANOTHER METHOD ***
        """
        assert self.contacts, "Need normal contacts first"
        
        self.logger.info("Mapping bbcontacts file to previously read contacts")
        # Read the contacts from the CASP RR formatted file. Unlike with other
        # contacts, bbcontacts contact pairs are also predicted around the turn
        # of a B-strand, so do not filter neighbours.
        bb_contacts = self._readContacts(secondContactfile, self.sequence)
                               
        self._iterBBcontactsIntoExistingContactList(bb_contacts)
                    
        return
    ##End _readAdditionalBBcontacts()
    
    def _iterBBcontactsIntoExistingContactList(self, bb_contacts):
        for bbcontact in bb_contacts:           
            isPresent=False
            for index, contact in enumerate(self.contacts):
                if bbcontact['res1_index'] == contact['res1_index'] and bbcontact['res2_index'] == contact['res2_index']:
                    self.contacts[index]['weight'] = 2
                    isPresent=True
                    
            if isPresent: self._checkBBcontactsNeighbours(bbcontact)
            else: self.contacts.append(bbcontact)
        return
    
    def _checkBBcontactsNeighbours(self, bbcontact):
        for index, contact in enumerate(self.contacts):
            for distance in [1, 2]:
                if contact['res1_index'] == (bbcontact['res1_index'] + distance) and contact['res2_index'] == bbcontact['res2_index']:
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == (bbcontact['res1_index'] - distance) and contact['res2_index'] == bbcontact['res2_index']:
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == bbcontact['res1_index'] and contact['res2_index'] == (bbcontact['res2_index']  + distance):
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == bbcontact['res1_index'] and contact['res2_index'] == (bbcontact['res2_index']  - distance):
                    self.contacts[index]['weight'] = 2 
        return
    
    def _readContacts(self, contactfile, sequence):
        '''Read the contactfile using the CASP RR Parser'''
        cp = parse_casprr.CaspContactParser()
        cp.read(contactfile)
        cp.sortContacts("raw_score", descending=True)
        cp.assignAminoAcids(sequence)
        cp.calculateScalarScores()
        return cp.contacts
    
    def _readConstraints(self, constraintsfile):
        '''Read the constraintsfile using parser'''
        cp = parse_constraints.ConstraintfileParser()
        try: 
            cp.read(constraintsfile)
        except ValueError, e: 
            self.logger.warning("Skipping reading of constraints:\n{0}".format(e))
            cp.contacts = []
        return cp.contacts
        

class Test(unittest.TestCase):
    def setUp(self):
        self.c = Contacter()
    
    def testFilter(self):
        contacts = [{'res1_index': 1, 'res2_index': 10},
                    {'res1_index': 1, 'res2_index': 6},
                    {'res1_index': 1, 'res2_index': 5},
                    {'res1_index': 1, 'res2_index': 2}]
        
        ref_contacts1 = [{'res1_index': 1, 'res2_index': 10},
                         {'res1_index': 1, 'res2_index': 6}]
        
        
        ref_contacts2 = [{'res1_index': 1, 'res2_index': 10},
                         {'res1_index': 1, 'res2_index': 6},
                         {'res1_index': 1, 'res2_index': 5},
                         {'res1_index': 1, 'res2_index': 2}]

        contacts_f1 = self.c._filterNeighbours(contacts, 5)
        contacts_f2 = self.c._filterNeighbours(contacts, 1)

        self.assertEqual(ref_contacts1, contacts_f1)
        self.assertEqual(ref_contacts2, contacts_f2)

    def testIter(self):
        self.c.contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                           {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                           {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                           {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                           {'res1_index': 5, 'res2_index':10, 'weight': 1}]
        
        bbcontacts = [{'res1_index':10, 'res2_index':100, 'weight': 1},
                      {'res1_index':11, 'res2_index': 99, 'weight': 1},
                      {'res1_index':12, 'res2_index': 98, 'weight': 1}]
        
        ref_contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                        {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                        {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                        {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                        {'res1_index': 5, 'res2_index':10, 'weight': 1},
                        {'res1_index':10, 'res2_index':100, 'weight': 1},
                        {'res1_index':11, 'res2_index':99, 'weight': 1},
                        {'res1_index':12, 'res2_index':98, 'weight': 1}]
        
        self.c._iterBBcontactsIntoExistingContactList(bbcontacts)
        
        self.assertItemsEqual(ref_contacts, self.c.contacts)        

    def testNeighbour(self):
        self.c.contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                           {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                           {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                           {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                           {'res1_index': 5, 'res2_index':10, 'weight': 1},
                    
                           {'res1_index': 3, 'res2_index': 6, 'weight': 1},
                           {'res1_index': 3, 'res2_index': 7, 'weight': 1},
                           {'res1_index': 3, 'res2_index': 9, 'weight': 1},
                           {'res1_index': 3, 'res2_index':10, 'weight': 1},
                    
                           {'res1_index': 1, 'res2_index': 8, 'weight': 1},
                           {'res1_index': 2, 'res2_index': 8, 'weight': 1},
                           {'res1_index': 4, 'res2_index': 8, 'weight': 1},
                           {'res1_index': 5, 'res2_index': 8, 'weight': 1},
                    
                           {'res1_index':10, 'res2_index':100, 'weight': 1}]
        
        bbcontact = {'res1_index': 3, 'res2_index': 8, 'weight': 1}
        
        ref_contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                        {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                        {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                        {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                        {'res1_index': 5, 'res2_index':10, 'weight': 1},
                    
                        {'res1_index': 3, 'res2_index': 6, 'weight': 2},
                        {'res1_index': 3, 'res2_index': 7, 'weight': 2},
                        {'res1_index': 3, 'res2_index': 9, 'weight': 2},
                        {'res1_index': 3, 'res2_index':10, 'weight': 2},
                    
                        {'res1_index': 1, 'res2_index': 8, 'weight': 2},
                        {'res1_index': 2, 'res2_index': 8, 'weight': 2},
                        {'res1_index': 4, 'res2_index': 8, 'weight': 2},
                        {'res1_index': 5, 'res2_index': 8, 'weight': 2},
                    
                        {'res1_index':10, 'res2_index':100, 'weight': 1}]
        
        
        self.c._checkBBcontactsNeighbours(bbcontact)
        
        self.assertItemsEqual(ref_contacts, self.c.contacts)

    def testPPV(self):
        contacts = [[19.089, 51.837, 57.937], [16.279, 47.215, 58.769],
                    [7.599, 48.942, 52.633], [14.616, 50.231, 60.890],
                    [6.783, 46.760, 42.296]]
        contacts = [numpy.array(x) for x in contacts]

        dist_mat = numpy.zeros((5, 5), numpy.float)
        dist_mat[0,0] = numpy.sqrt(numpy.sum((contacts[0] - contacts[0]) ** 2))
        dist_mat[0,1] = numpy.sqrt(numpy.sum((contacts[0] - contacts[1]) ** 2))
        dist_mat[0,2] = numpy.sqrt(numpy.sum((contacts[0] - contacts[2]) ** 2))
        dist_mat[0,3] = numpy.sqrt(numpy.sum((contacts[0] - contacts[3]) ** 2))
        dist_mat[0,4] = numpy.sqrt(numpy.sum((contacts[0] - contacts[4]) ** 2))
        dist_mat[1,0] = numpy.sqrt(numpy.sum((contacts[1] - contacts[0]) ** 2))
        dist_mat[1,1] = numpy.sqrt(numpy.sum((contacts[1] - contacts[1]) ** 2))
        dist_mat[1,2] = numpy.sqrt(numpy.sum((contacts[1] - contacts[2]) ** 2))
        dist_mat[1,3] = numpy.sqrt(numpy.sum((contacts[1] - contacts[3]) ** 2))
        dist_mat[1,4] = numpy.sqrt(numpy.sum((contacts[1] - contacts[4]) ** 2))
        dist_mat[2,0] = numpy.sqrt(numpy.sum((contacts[2] - contacts[0]) ** 2))
        dist_mat[2,1] = numpy.sqrt(numpy.sum((contacts[2] - contacts[1]) ** 2))
        dist_mat[2,2] = numpy.sqrt(numpy.sum((contacts[2] - contacts[2]) ** 2))
        dist_mat[2,3] = numpy.sqrt(numpy.sum((contacts[2] - contacts[3]) ** 2))
        dist_mat[2,4] = numpy.sqrt(numpy.sum((contacts[2] - contacts[4]) ** 2))
        dist_mat[3,0] = numpy.sqrt(numpy.sum((contacts[3] - contacts[0]) ** 2))
        dist_mat[3,1] = numpy.sqrt(numpy.sum((contacts[3] - contacts[1]) ** 2))
        dist_mat[3,2] = numpy.sqrt(numpy.sum((contacts[3] - contacts[2]) ** 2))
        dist_mat[3,3] = numpy.sqrt(numpy.sum((contacts[3] - contacts[3]) ** 2))
        dist_mat[3,4] = numpy.sqrt(numpy.sum((contacts[3] - contacts[4]) ** 2))
        dist_mat[4,0] = numpy.sqrt(numpy.sum((contacts[4] - contacts[0]) ** 2))
        dist_mat[4,1] = numpy.sqrt(numpy.sum((contacts[4] - contacts[1]) ** 2))
        dist_mat[4,2] = numpy.sqrt(numpy.sum((contacts[4] - contacts[2]) ** 2))
        dist_mat[4,3] = numpy.sqrt(numpy.sum((contacts[4] - contacts[3]) ** 2))
        dist_mat[4,4] = numpy.sqrt(numpy.sum((contacts[4] - contacts[4]) ** 2))

        RCmap = dist_mat < 8
        RCmap_seq = "AAAAA"

        input_contacts1 = [{'res1_index': 1, 'res2_index': 3},
                           {'res1_index': 1, 'res2_index': 4},
                           {'res1_index': 2, 'res2_index': 4},
                           {'res1_index': 3, 'res2_index': 4}]

        input_contacts2 = [{'res1_index': 1, 'res2_index': 2},
                           {'res1_index': 1, 'res2_index': 4},
                           {'res1_index': 2, 'res2_index': 4},
                           {'res1_index': 5, 'res2_index': 5}]

        input_contacts3 = [{'res1_index': 1, 'res2_index': 3},
                           {'res1_index': 2, 'res2_index': 5},
                           {'res1_index': 3, 'res2_index': 4}]

        input_contacts4 = [{'res1_index': 1, 'res2_index': 2},
                           {'res1_index': 1, 'res2_index': 3},
                           {'res1_index': 2, 'res2_index': 4},
                           {'res1_index': 2, 'res2_index': 5},
                           {'res1_index': 3, 'res2_index': 4}]

        ppv1 = self.c._ppv_score(input_contacts1, RCmap, RCmap_seq)
        ppv2 = self.c._ppv_score(input_contacts2, RCmap, RCmap_seq)
        ppv3 = self.c._ppv_score(input_contacts3, RCmap, RCmap_seq)
        ppv4 = self.c._ppv_score(input_contacts4, RCmap, RCmap_seq)

        self.assertEqual(0.5, ppv1)
        self.assertEqual(1.0, ppv2)
        self.assertEqual(0.0, ppv3)
        self.assertEqual(0.4, ppv4)

    def testRefContacts(self):
        contacts = [[19.089, 51.837, 57.937], [16.279, 47.215, 58.769],
                    [7.599, 48.942, 52.633], [14.616, 50.231, 60.890],
                    [6.783, 46.760, 42.296]]
        contacts = [numpy.array(x) for x in contacts]

        length = 5
        dist_mat = numpy.zeros((length, length), numpy.float)
        dist_mat[0,0] = numpy.sqrt(numpy.sum((contacts[0] - contacts[0]) ** 2))
        dist_mat[0,1] = numpy.sqrt(numpy.sum((contacts[0] - contacts[1]) ** 2))
        dist_mat[0,2] = numpy.sqrt(numpy.sum((contacts[0] - contacts[2]) ** 2))
        dist_mat[0,3] = numpy.sqrt(numpy.sum((contacts[0] - contacts[3]) ** 2))
        dist_mat[0,4] = numpy.sqrt(numpy.sum((contacts[0] - contacts[4]) ** 2))
        dist_mat[1,0] = numpy.sqrt(numpy.sum((contacts[1] - contacts[0]) ** 2))
        dist_mat[1,1] = numpy.sqrt(numpy.sum((contacts[1] - contacts[1]) ** 2))
        dist_mat[1,2] = numpy.sqrt(numpy.sum((contacts[1] - contacts[2]) ** 2))
        dist_mat[1,3] = numpy.sqrt(numpy.sum((contacts[1] - contacts[3]) ** 2))
        dist_mat[1,4] = numpy.sqrt(numpy.sum((contacts[1] - contacts[4]) ** 2))
        dist_mat[2,0] = numpy.sqrt(numpy.sum((contacts[2] - contacts[0]) ** 2))
        dist_mat[2,1] = numpy.sqrt(numpy.sum((contacts[2] - contacts[1]) ** 2))
        dist_mat[2,2] = numpy.sqrt(numpy.sum((contacts[2] - contacts[2]) ** 2))
        dist_mat[2,3] = numpy.sqrt(numpy.sum((contacts[2] - contacts[3]) ** 2))
        dist_mat[2,4] = numpy.sqrt(numpy.sum((contacts[2] - contacts[4]) ** 2))
        dist_mat[3,0] = numpy.sqrt(numpy.sum((contacts[3] - contacts[0]) ** 2))
        dist_mat[3,1] = numpy.sqrt(numpy.sum((contacts[3] - contacts[1]) ** 2))
        dist_mat[3,2] = numpy.sqrt(numpy.sum((contacts[3] - contacts[2]) ** 2))
        dist_mat[3,3] = numpy.sqrt(numpy.sum((contacts[3] - contacts[3]) ** 2))
        dist_mat[3,4] = numpy.sqrt(numpy.sum((contacts[3] - contacts[4]) ** 2))
        dist_mat[4,0] = numpy.sqrt(numpy.sum((contacts[4] - contacts[0]) ** 2))
        dist_mat[4,1] = numpy.sqrt(numpy.sum((contacts[4] - contacts[1]) ** 2))
        dist_mat[4,2] = numpy.sqrt(numpy.sum((contacts[4] - contacts[2]) ** 2))
        dist_mat[4,3] = numpy.sqrt(numpy.sum((contacts[4] - contacts[3]) ** 2))
        dist_mat[4,4] = numpy.sqrt(numpy.sum((contacts[4] - contacts[4]) ** 2))

        ref_map = dist_mat < 8
        ref_contacts = numpy.where(dist_mat < 8)
        
        contacts_out, contact_map = self.c._cb_contacts(contacts, contacts, length)
        
        numpy.testing.assert_equal(ref_contacts, contacts_out)
        numpy.testing.assert_equal(ref_map, contact_map)

    def testTruncate(self):
        contacts = [{'res1_index': 1, 'res2_index': 10},
                    {'res1_index': 1, 'res2_index': 6},
                    {'res1_index': 1, 'res2_index': 5},
                    {'res1_index': 1, 'res2_index': 2}]
        
        ref_c1 = [{'res1_index': 1, 'res2_index': 10},
                  {'res1_index': 1, 'res2_index': 6},
                  {'res1_index': 1, 'res2_index': 5},
                  {'res1_index': 1, 'res2_index': 2}]
        ref_c2 = [{'res1_index': 1, 'res2_index': 10},
                  {'res1_index': 1, 'res2_index': 6}]
        ref_c3 = []
        ref_c4 = [{'res1_index': 1, 'res2_index': 10},
                  {'res1_index': 1, 'res2_index': 6},
                  {'res1_index': 1, 'res2_index': 5},
                  {'res1_index': 1, 'res2_index': 2}]
        
        out_c1 = self.c._truncateContactList(contacts, int(4*1.0))
        out_c2 = self.c._truncateContactList(contacts, int(4*0.7))
        out_c3 = self.c._truncateContactList(contacts, int(4*0.0))
        out_c4 = self.c._truncateContactList(contacts, int(4*2.0))        
        
        self.assertEqual(ref_c1, out_c1)
        self.assertEqual(ref_c2, out_c2)
        self.assertEqual(ref_c3, out_c3)
        self.assertEqual(ref_c4, out_c4)

    

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=str, default=None, dest='bbcontacts_file',
                        help="Additional bbcontacts CASPRR contactfile")
    parser.add_argument('-c', type=float, default=1.0, dest="constraints_factor",
                        help="Defines number of contacts to use (L/*x*)")
    parser.add_argument('-d', type=int, default=5, dest="distance_to_neighbour",
                        help="Defines sequence cutoff to neighbouring residues")
    parser.add_argument('-e', type=str, default="FADE", dest='energy_function',
                        help="Rosetta function to use")
    parser.add_argument('fasta')
    parser.add_argument('-n', type=str, default="ampl_", dest="name",
                        help="Job name")
    parser.add_argument('-o', type=str, default="ampl_.cst", dest="output",
                        help="Output file")
    parser.add_argument('-s', type=str, default=None, dest="native_pdb",
                        help="Reference structure")
    parser.add_argument('-ss2', type=str, default=None, dest="psipred_ss2",
                        help="Secondary structure prediction")
    
    contacts = parser.add_mutually_exclusive_group(required=True)
    contacts.add_argument('-contact_file', default=None)
    contacts.add_argument('-constraints_file', default=None)
    
    option = parser.add_mutually_exclusive_group(required=True)
    option.add_argument('-format', action="store_true")
    option.add_argument('-plot', action="store_true")
    option.add_argument('-ppv', action="store_true")

    optd = vars(parser.parse_args())
    
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    # Set working directory
    optd['work_dir'] = os.getcwd()

    # Reformat to what we need
    logging.debug('Parsing FASTA file')
    try: fp = ample_sequence.Sequence(fasta=optd['fasta'])
    except Exception as e:
        print "Error parsing FASTA file: {0}\n\n{1}".format(optd['fasta'],e.message)
    if fp.numSequences() != 1:
        print "ERROR! Fasta file {0} has > 1 sequence in it.".format(optd['fasta'])
    optd['sequence'] = fp.sequence()

    c = Contacter(optd)
    if optd['format']: c.format(optd['output'])
    elif optd['plot']: c.plot("ampl_.cm.pdf", ss2file=optd['psipred_ss2'], structurefile=optd['native_pdb'])
    elif optd['ppv']: c.ppv(optd['native_pdb'])
    else: logging.critical("No option selected")
