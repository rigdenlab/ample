# # #                                                # # #
# This file contains Rosetta energy function definitions #
# for constraint usage. Different functions are defined  #
# and additional ones can be added with set parameters.  #
#                                                        #
# All functions take a formatted contact dictionary as   #
# input. Ensure parameters match the keys in the dict.   #
# # #                                                # # #

def FADE(contact):
    template="AtomPair %(atom1)s %(res1_index)d %(atom2)s %(res2_index)d FADE -10 19 10 %(weight).2f 0" 
    
    return template % {'atom1': contact['atom1'], 'res1_index': contact['res1_index'],
                       'atom2': contact['atom2'], 'res2_index': contact['res2_index'],
                       'weight': (contact['weight'] * -15.00)}


def FADE_default(contact):  
    template = "AtomPair %(atom1)s %(res1_index)s %(atom2)s %(res2_index)s FADE -10 19 10 -15.00 0"
    
    return template % {'atom1': contact['atom1'], 'res1_index': contact['res1_index'],
                       'atom2': contact['atom2'], 'res2_index': contact['res2_index']}

