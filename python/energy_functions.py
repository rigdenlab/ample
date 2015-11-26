# # #                                                # # #
# This file contains Rosetta energy function definitions #
# for constraint usage. Different functions are defined  #
# and additional ones can be added with set parameters.  #
#                                                        #
# All functions take a formatted contact dictionary as   #
# input. Ensure parameters match the keys in the dict.   #
# # #                                                # # #

# # #
#
# A contact input dictionary contains the following keys.
# Code defining a contact in parse_contactfile.         
#
# MAIN KEYS:
#
# `atom1`            -> Stores the ATOM of residue 1
# `atom2`            -> Stores the ATOM of residue 2
# `res1`             -> Stores the AMINO ACID of residue 1
# `res2`             -> Stores the AMINO ACID of residue 2
# `res1_index`       -> Stores the RESIDUE INDEX of residue1
# `res2_index`       -> Stores the RESIDUE INDEX of residue2
# `confidence_score` -> Stores the CONTACT PREDICTOR confidence value
# `weight`           -> Stores a WEIGHT factor, i.e number of times a contact was predicted
# `lb`               -> Stores the LB value of a contact, CASP RR 3rd column
# `ub`               -> Stores the UB value of a contact, CASP RR 4th column
# `true_positive`    -> Stores TRUE by default and only changed if a native structure is provided
# `method`           -> Stores the METHOD origin of the contacts
# `file`             -> Stores the FILE where contacts were stored
#
# # #
#
# BBCONTACTS-SPECIFIC KEYS:
#
# `diversity_factor`        -> Stores the diversity factor used
# `internal_strand_position -> Stores the internal strand position
# `strand_index`            -> Stores the index of the strand
# `strand_orientation`      -> Stores the orientation of the a strand
#
# # #


def FADE(contact):
    template="AtomPair %(atom1)s %(res1_index)d %(atom2)s %(res2_index)d FADE -10 19 10 %(weight).2f 0" 
    
    return template % {'atom1': contact['atom1'], 'res1_index': contact['res1_index'],
                       'atom2': contact['atom2'], 'res2_index': contact['res2_index'],
                       'weight': (contact['weight'] * -15.00)}


def FADE_default(contact):  
    template = "AtomPair %(atom1)s %(res1_index)s %(atom2)s %(res2_index)s FADE -10 19 10 -15.00 0"
    
    return template % {'atom1': contact['atom1'], 'res1_index': contact['res1_index'],
                       'atom2': contact['atom2'], 'res2_index': contact['res2_index']}

