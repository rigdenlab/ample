"""Test functions for python.contacts"""

import numpy
import unittest

from ample.python.contacts import Contacter

class Test(unittest.TestCase):
    def setUp(self):
        self.c = Contacter()
    
    def test_filter(self):
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

    def test_iter(self):
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

    def test_neighbour(self):
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

    def test_PPV(self):
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

    def test_refContacts(self):
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

    def test_truncate(self):
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
    unittest.main()    
