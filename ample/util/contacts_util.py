#!/usr/bin/env ccp4-python

import logging
import numpy
import os

from ample.modelling import energy_functions
from ample.parsers import casprr_parser
from ample.parsers import psipred_parser
from ample.parsers import restraints_parser
from ample.util import pdb_edit

__author__ = "Felix Simkovic"
__date__ = "14.11.2015"
__version__ = "1.0"

# Check matplotlib is available - problems with Windows
try:
    from ample.util import plot_util
    _MATPLOTLIB = True
except (ImportError, ValueError):
    _MATPLOTLIB = False

# Check benchmark_util is available to standardise PDB if necessary
try:
    from ample.util import benchmark_util
    _BENCHMARK = True
except ImportError:
    _BENCHMARK = False

# Check biopython is availbale - compatibility with older ccp4-python 
try:
    from ample.parsers import alignment_parser
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False

LOGGER = logging.getLogger(__name__)

def checkOptions(amoptd):
    """Function to check that all contact files are available"""
    
    # Make sure contact file is provided with bbcontacts_file
    if not amoptd['contact_file'] and amoptd['bbcontacts_file']:
        msg = "Must provide -contact_file when using -bbcontacts_file or use as -contact_file"
        LOGGER.critical(msg)
        raise RuntimeError(msg)
    
    # Check the existence of the contact file and whether it is CASP RR format
    # based on the required header `PFRMAT RR`
    if amoptd['contact_file']:
        if not os.path.exists(str(amoptd['contact_file'])):
            msg = "Cannot find contact file:\n{0}".format(amoptd['contact_file'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
           
        if not casprr_parser.CaspContactParser().check_format(amoptd['contact_file']):
            msg = "Wrong format in contact file:\n{0}".format(amoptd['contact_file'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
    
    # Check the existence of the contact file and whether it is CASP RR format
    # based on the required header `PFRMAT RR`   
    if amoptd['bbcontacts_file']:
        if not os.path.exists(amoptd['bbcontacts_file']):
            msg = "Cannot find contact file:\n{0}".format(amoptd['contact_file'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
            
        if not casprr_parser.CaspContactParser().check_format(amoptd['bbcontacts_file']):
            msg = "Wrong format in contact file:\n{0}".format(amoptd['bbcontacts_file'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
    
    # Make sure user selected energy function is pre-defined
    if amoptd['energy_function']:
        try: 
            # We can discard the returned value as we only want to check if it
            # exists in the module
            _ = getattr(energy_functions, amoptd['energy_function'])
        except AttributeError:
            msg = "Rosetta energy function {0} unavailable".format(amoptd['energy_function'])
            LOGGER.critical(msg)
            raise RuntimeError(msg)
    return


class Contacter(object):
    """Class to handle contact predictions"""
    
    def __init__(self, optd=None):
        
        self.bbcontacts_file = None
        self.contact_file = None
        self.contact_map = None
        self.contact_ppv = None
        self.distance_to_neighbour = 5
        self.native_cutoff = 8
        self.psipred_ss2 = None
        self.restraints_file = None
        self.restraints_factor = 1.0
        self.restraints_infile = None
        self.structure_pdb = None
        
        self.contacts = None
        self.raw_contacts = None
        self.sequence = None
        
        self.energy_function = None
        
        if optd: self.init(optd)
        
        return
    
    def init(self, optd):

        self.contact_file = optd['contact_file']
        self.contact_map = os.path.join(optd['work_dir'], optd['name'] + ".cm.pdf")
        self.restraints_file = os.path.join(optd['work_dir'], optd['name'] + ".cst")
        self.restraints_infile = optd['restraints_file']
        self.sequence = optd['sequence']
        
        self.energy_function = optd['energy_function']
        
        # Optional files
        if optd['native_pdb'] and optd['native_pdb_std']: 
            self.structure_pdb = optd['native_pdb_std']
        elif optd['native_pdb'] and _BENCHMARK:
            benchmark_util.analysePdb(optd)
            self.structure_pdb = optd['native_pdb_std']
        if optd['native_cutoff']: 
            self.native_cutoff = optd['native_cutoff']
        if optd['psipred_ss2']: 
            self.psipred_ss2 = optd['psipred_ss2']
        
        # Set these parameters class wide so we can use them in the plot function
        # in case a second map is given
        self.distance_to_neighbour = optd['distance_to_neighbour']
        self.restraints_factor = optd['restraints_factor']
        
        # Extract the raw sequence
        if self.contact_file:
            self.raw_contacts = self._readContacts(self.contact_file, self.sequence)
            # Prepare the contacts depending on all the user/default options
            self.contacts = self._prepare(self.raw_contacts, 
                                          self.restraints_factor, 
                                          self.distance_to_neighbour)
            
            # Map bbcontacts on top of the previously provided contacts
            if optd['bbcontacts_file']: 
                self.bbcontacts_file = optd['bbcontacts_file']
                self._readAdditionalBBcontacts(optd['bbcontacts_file'])
                
        elif self.restraints_infile:
            self.contacts = self._readRestraints(self.restraints_infile)
            
        return
                
    def format(self, restraintfile):
        """Format contacts to Rosetta restraints"""
        
        LOGGER.info("Re-formatting contacts to restraints " +
                    "using the {0} function".format(self.energy_function))
        
        # Format the contacts to restraints
        contact_formatted_lines = self._formatToRestraints(self.contacts, 
                                                           self.energy_function)
        
        # Write contacts to restraint file
        with open(restraintfile, 'w') as oh: 
            oh.write("\n".join(contact_formatted_lines))
        
        return
       
    def _formatToRestraints(self, contacts, user_function):
        """Return a list of Rosetta string lines"""
    
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
        self.format(self.restraints_file)
           
        # Contact map plotting. We can parse ss2file and structure file blindly, checks in place
        self.plot(self.contact_map,
                  ss2file=self.psipred_ss2,
                  structurefile=self.structure_pdb)
        
        if self.structure_pdb: 
            self.contact_ppv=self.ppv(self.structure_pdb)
        
        return
    
    def process_restraintsfile(self):
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
        ## MOVE ME - fix for GREMLIN prediction with commented bits on top
        with open(self.restraints_infile, "r") as in_fh, \
             open(self.restraints_file, "w") as out_fh:
            for line in iter(in_fh.readline, ''):
                if not line.strip() or line.startswith("#"): continue
                out_fh.write(line)

        return
          
    def plot(self, figurefile, second_contact_file=None, ss2file=None, structurefile=None, offset=0):
        """Plot a contact map"""
        
        # make sure we can import matplotlib
        if not _MATPLOTLIB: 
            LOGGER.warning("Cannot plot contact map due to missing Python dependencies: MATPLOTLIB")
            return
        
        # Just to make sure we have either a structurefile or a pre-defined structure
        # and we can import Biopython
        if structurefile and _BIOPYTHON:
            available_structure = True
        elif self.structure_pdb and _BIOPYTHON:
            available_structure = True
            structurefile = self.structure_pdb
        else:
            available_structure = False
            
        ap = plot_util.Plotter()
        ap.initialise(figsize=(3,3), dpi=600)
        
        # Get the coordinates from the reference structure and plot in gray
        if available_structure:
            RCs, RCm = self.structureContacts(structurefile, self.sequence)
            RCs = list([i+offset for i in R] for R in RCs)
            ap.addScatter(RCs[0], RCs[1], marker='.', c='#DDDDDD', s=10,
                          edgecolor="none", linewidths=0.0)
        else:
            RCs, RCm = numpy.empty(shape=(0,0)), numpy.empty(shape=(0,0))

        # Get the coordinates for the secondary structure prediction and plot it alogn the diagonal
        if ss2file:
            SCs, Scolors = self.secondaryStructureContacts(ss2file)
            SCs = list([i+offset for i in S] for S in SCs)
            ap.addScatter(SCs[0], SCs[1], marker='.', c=Scolors, s=35,
                          edgecolor="black", linewidths=0.1)

        # Bit cleaner code if we extract X and Y coordinates and then parse them two the plotter
        # Reduce residue index by one to account for counting from 0
        Xs = [ i['res1_index']+offset-1 for i in self.contacts ]
        Ys = [ i['res2_index']+offset-1 for i in self.contacts ]
        
        # Get the color codes for main contact map
        color_codes = self.get_color_codes(self.contacts, RCm=RCm, offset=offset)
        
        # Allow the plotting of two contact maps if a second set of contacts
        # is provided
        if second_contact_file:
            second_raw_contacts = self._readContacts(second_contact_file,
                                                     self.sequence)
            # Prepare the contacts depending on all the user/default options
            second_contacts = self._prepare(second_raw_contacts, 
                                            self.restraints_factor, 
                                            self.distance_to_neighbour)
            color_codes_2 = self.get_color_codes(second_contacts, RCm=RCm, offset=offset)
            X2s = [ i['res1_index']+offset-1 for i in second_contacts ]
            Y2s = [ i['res2_index']+offset-1 for i in second_contacts ]
        else:
            X2s, Y2s = Xs, Ys
            color_codes_2 = color_codes
        
        ap.addScatter(Xs, Ys, marker='.', c=color_codes, s=10,
                      edgecolor="none", linewidths=0.0)
        ap.addScatter(Y2s, X2s, marker='.', c=color_codes_2, s=10,
                      edgecolor="none", linewidths=0.0)
        
        ap.axisLimits(len(self.sequence), offset=offset)
        #ap.axisTitles("Residue number", "Residue number")

        ap.saveFig(figurefile)
        
        return
    
    def get_color_codes(self, contacts, RCm=None, offset=0):
        """Wrapper for getting the color codes"""
        return self._color_codes(contacts, RCm=RCm, offset=offset)
    
    def _color_codes(self, contacts, RCm, offset):
        '''Get the color codes for each contact depending on match and weight'''
        
        # Definitions of colors
        color_choices = {'true_positive_double': '#2D9D00',
                         'true_positive_single': '#38C700',
                         'false_positive_double': '#AB0000',
                         'false_positive_single': '#D70909',
                         'default_double': '#004F9D',
                         'default_single': '#0482ff',
        }
        
        # Store the color codes per residue contact
        contact_colors = []
        
        for idx in range(len(contacts)):
            c_x = contacts[idx]['res1_index']-offset-1
            c_y = contacts[idx]['res2_index']-offset-1

            # No structure, use default colors but respect single/double weight
            if not RCm.any():
                color = color_choices['default_double'] \
                            if contacts[idx]['weight'] == 2 else \
                                color_choices['default_single']
            
            # Determine True Positives based on structure
            elif RCm[c_x, c_y] > 0:     
                color = color_choices['true_positive_double'] \
                            if contacts[idx]['weight'] == 2 else \
                                color_choices['true_positive_single']
                                
            # Determine False Positives based on structure
            elif RCm[c_x, c_y] == 0:
                color = color_choices['false_positive_double'] \
                            if contacts[idx]['weight'] == 2 else \
                                color_choices['false_positive_single']
                                
            # Make sure we have a color assigned
            if not color:
                msg = "Cannot handle a contact: {0} - {1}".format(c_x, c_y)
                raise Exception(msg)
            
            # Save the color for this particular contact
            contact_colors.append(color)
                
        return contact_colors
    
    def ppv(self, structurefile):
        """Calculate the True Positive Rate of contact prediction"""
         
        # Check that we were able to import Biopython
        if not _BIOPYTHON: return 0.0
        
        assert structurefile, "You need to provide PDB to calculate the PPV"
        
        _, RCm = self.structureContacts(structurefile, self.sequence)
        
        ppv = self._ppv_score(self.contacts, RCm, self.structure_seq)
        
        LOGGER.info("Accuracy of contact prediction (PPV): {0} %".format(ppv*100))
        
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
        pp = psipred_parser.PsipredSs2Parser(ss2file)
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
            aligned_seq_list = alignment_parser.AlignmentParser().align_sequences(alignmentSequence, 
                                                                                  self.structure_seq)
   
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
        
        contacts, contact_map = self._cb_contacts(cb_lst, cb_lst, len(self.structure_seq), self.native_cutoff)

        return contacts, contact_map
    
    def _cb_contacts(self, cb1_lst, cb2_lst, length, cutoff=8):
        '''Get the contacts between the two lists of contacts'''
        
        LOGGER.info("Distance cutoff of participating atoms is: %.1fA" % cutoff)

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
    
    def _prepare(self, contacts, restraint_factor, distance_to_neighbour):
        # No processing so far, but we will need to do it to match atom to contact
        nrRestraints = int(len(self.sequence) * restraint_factor)
        
        # Filter out any contact pairs that are within `distance residues` of each other
        contacts_filtered = self._filterNeighbours(contacts, distance_to_neighbour)
        
        # Truncate to the number of restraints defined by the user
        contacts_filtered_truncated = self._truncateContactList(contacts_filtered, nrRestraints)
        
        return contacts_filtered_truncated
   
    def _filterNeighbours(self, contacts, distance):
        ''' Filter any contact pairs closer than `distance` residues '''
        contacts_filtered = [contact for contact in contacts \
                                if abs(contact['res2_index'] - contact['res1_index']) >= distance]
        return contacts_filtered  
    
    def _truncateContactList(self, contacts, nrRestraints):
        return contacts[:nrRestraints]
    
    def _readAdditionalBBcontacts(self, secondContactfile):
        """ Method designated solely for reading bbcontacts to overlay/match
            them to the existing set of contacts.
            
            *** BE CAREFUL WHEN USING WITH ANOTHER METHOD ***
        """
        assert self.contacts, "Need normal contacts first"
        
        LOGGER.info("Mapping bbcontacts file to previously read contacts")
        # Read the contacts from the CASP RR formatted file. Unlike with other
        # contacts, bbcontacts contact pairs are also predicted around the turn
        # of a B-strand, so do not filter neighbours.
        bb_contacts = self._readContacts(secondContactfile, self.sequence)
                               
        self._iterBBcontactsIntoExistingContactList(bb_contacts)
                    
        return
    
    def _iterBBcontactsIntoExistingContactList(self, bb_contacts):
        for bbcontact in bb_contacts:           
            isPresent=False
            for index, contact in enumerate(self.contacts):
                if bbcontact['res1_index'] == contact['res1_index'] and \
                   bbcontact['res2_index'] == contact['res2_index']:
                    self.contacts[index]['weight'] = 2
                    isPresent=True
                    
            if isPresent: self._checkBBcontactsNeighbours(bbcontact)
            else: self.contacts.append(bbcontact)
        return
    
    def _checkBBcontactsNeighbours(self, bbcontact):
        for index, contact in enumerate(self.contacts):
            for distance in [1, 2]:
                if contact['res1_index'] == (bbcontact['res1_index'] + distance) and \
                   contact['res2_index'] == bbcontact['res2_index']:
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == (bbcontact['res1_index'] - distance) and \
                   contact['res2_index'] == bbcontact['res2_index']:
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == bbcontact['res1_index'] and \
                   contact['res2_index'] == (bbcontact['res2_index']  + distance):
                    self.contacts[index]['weight'] = 2 
                if contact['res1_index'] == bbcontact['res1_index'] and \
                   contact['res2_index'] == (bbcontact['res2_index']  - distance):
                    self.contacts[index]['weight'] = 2 
        return
    
    def _readContacts(self, contactfile, sequence):
        '''Read the contactfile using the CASP RR Parser'''
        cp = casprr_parser.CaspContactParser()
        cp.read(contactfile)
        cp.sort_contacts("raw_score", descending=True)
        cp.assign_amino_acids(sequence)
        cp.calculate_scalar_scores()
        return cp.contacts
    
    def _readRestraints(self, restraintsfile):
        '''Read the restraintsfile using parser'''
        cp = restraints_parser.RestraintfileParser()
        try: 
            cp.read(restraintsfile)
        except ValueError, e: 
            LOGGER.warning("Skipping reading of restraints for processing:\n{0}".format(e))
            cp.contacts = []
        return cp.contacts
        
