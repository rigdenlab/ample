"""Test functions for util.contacts_util"""

import numpy
import unittest

from ample.util import contacts_util

class Test(unittest.TestCase):
    
    def test_filter(self):
        input = [{'res1_index': 1, 'res2_index': 10},
                 {'res1_index': 1, 'res2_index': 6},
                 {'res1_index': 1, 'res2_index': 5},
                 {'res1_index': 1, 'res2_index': 2}]
        c = contacts_util.Contacter()
        input_f = c._filterNeighbours(input, 5)
        ref_output = [{'res1_index': 1, 'res2_index': 10},
                      {'res1_index': 1, 'res2_index': 6}]
        self.assertEqual(ref_output, input_f)
        input_f = c._filterNeighbours(input, 1)
        ref_output = [{'res1_index': 1, 'res2_index': 10},
                        {'res1_index': 1, 'res2_index': 6},
                        {'res1_index': 1, 'res2_index': 5},
                        {'res1_index': 1, 'res2_index': 2}]
        self.assertEqual(ref_output, input_f)

    def test_iter(self):
        c = contacts_util.Contacter()
        c.contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                      {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                      {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                      {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                      {'res1_index': 5, 'res2_index':10, 'weight': 1}]
        bbcontacts = [{'res1_index':10, 'res2_index':100, 'weight': 1},
                      {'res1_index':11, 'res2_index': 99, 'weight': 1},
                      {'res1_index':12, 'res2_index': 98, 'weight': 1}]
        c._iterBBcontactsIntoExistingContactList(bbcontacts)
        ref_contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
                        {'res1_index': 2, 'res2_index': 7, 'weight': 1},
                        {'res1_index': 3, 'res2_index': 8, 'weight': 1},
                        {'res1_index': 4, 'res2_index': 9, 'weight': 1},
                        {'res1_index': 5, 'res2_index':10, 'weight': 1},
                        {'res1_index':10, 'res2_index':100, 'weight': 1},
                        {'res1_index':11, 'res2_index':99, 'weight': 1},
                        {'res1_index':12, 'res2_index':98, 'weight': 1}]
        self.assertItemsEqual(ref_contacts, c.contacts)        

    def test_neighbour(self):
        c = contacts_util.Contacter()
        c.contacts = [{'res1_index': 1, 'res2_index': 6, 'weight': 1},
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
        c._checkBBcontactsNeighbours(bbcontact)
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
        self.assertItemsEqual(ref_contacts, c.contacts)

    def test_PPV(self):
        input = [[19.089, 51.837, 57.937], [16.279, 47.215, 58.769],
                 [7.599, 48.942, 52.633], [14.616, 50.231, 60.890],
                 [6.783, 46.760, 42.296]]
        input = [numpy.array(x) for x in input]

        dist_mat = numpy.zeros((5, 5), numpy.float)
        dist_mat[0,0] = numpy.sqrt(numpy.sum((input[0] - input[0]) ** 2))
        dist_mat[0,1] = numpy.sqrt(numpy.sum((input[0] - input[1]) ** 2))
        dist_mat[0,2] = numpy.sqrt(numpy.sum((input[0] - input[2]) ** 2))
        dist_mat[0,3] = numpy.sqrt(numpy.sum((input[0] - input[3]) ** 2))
        dist_mat[0,4] = numpy.sqrt(numpy.sum((input[0] - input[4]) ** 2))
        dist_mat[1,0] = numpy.sqrt(numpy.sum((input[1] - input[0]) ** 2))
        dist_mat[1,1] = numpy.sqrt(numpy.sum((input[1] - input[1]) ** 2))
        dist_mat[1,2] = numpy.sqrt(numpy.sum((input[1] - input[2]) ** 2))
        dist_mat[1,3] = numpy.sqrt(numpy.sum((input[1] - input[3]) ** 2))
        dist_mat[1,4] = numpy.sqrt(numpy.sum((input[1] - input[4]) ** 2))
        dist_mat[2,0] = numpy.sqrt(numpy.sum((input[2] - input[0]) ** 2))
        dist_mat[2,1] = numpy.sqrt(numpy.sum((input[2] - input[1]) ** 2))
        dist_mat[2,2] = numpy.sqrt(numpy.sum((input[2] - input[2]) ** 2))
        dist_mat[2,3] = numpy.sqrt(numpy.sum((input[2] - input[3]) ** 2))
        dist_mat[2,4] = numpy.sqrt(numpy.sum((input[2] - input[4]) ** 2))
        dist_mat[3,0] = numpy.sqrt(numpy.sum((input[3] - input[0]) ** 2))
        dist_mat[3,1] = numpy.sqrt(numpy.sum((input[3] - input[1]) ** 2))
        dist_mat[3,2] = numpy.sqrt(numpy.sum((input[3] - input[2]) ** 2))
        dist_mat[3,3] = numpy.sqrt(numpy.sum((input[3] - input[3]) ** 2))
        dist_mat[3,4] = numpy.sqrt(numpy.sum((input[3] - input[4]) ** 2))
        dist_mat[4,0] = numpy.sqrt(numpy.sum((input[4] - input[0]) ** 2))
        dist_mat[4,1] = numpy.sqrt(numpy.sum((input[4] - input[1]) ** 2))
        dist_mat[4,2] = numpy.sqrt(numpy.sum((input[4] - input[2]) ** 2))
        dist_mat[4,3] = numpy.sqrt(numpy.sum((input[4] - input[3]) ** 2))
        dist_mat[4,4] = numpy.sqrt(numpy.sum((input[4] - input[4]) ** 2))

        RCmap = dist_mat < 8
        RCmap_seq = "AAAAA"

        c = contacts_util.Contacter()
        input_contacts = [{'res1_index': 1, 'res2_index': 3},
                          {'res1_index': 1, 'res2_index': 4},
                          {'res1_index': 2, 'res2_index': 4},
                          {'res1_index': 3, 'res2_index': 4}]
        ppv = c._ppv_score(input_contacts, RCmap, RCmap_seq)
        self.assertEqual(0.5, ppv)
        
        input_contacts = [{'res1_index': 1, 'res2_index': 2},
                          {'res1_index': 1, 'res2_index': 4},
                          {'res1_index': 2, 'res2_index': 4},
                          {'res1_index': 5, 'res2_index': 5}]
        ppv = c._ppv_score(input_contacts, RCmap, RCmap_seq)
        self.assertEqual(1.0, ppv)
        
        input_contacts = [{'res1_index': 1, 'res2_index': 3},
                          {'res1_index': 2, 'res2_index': 5},
                          {'res1_index': 3, 'res2_index': 4}]
        ppv = c._ppv_score(input_contacts, RCmap, RCmap_seq)
        self.assertEqual(0.0, ppv)
        
        input_contacts = [{'res1_index': 1, 'res2_index': 2},
                          {'res1_index': 1, 'res2_index': 3},
                          {'res1_index': 2, 'res2_index': 4},
                          {'res1_index': 2, 'res2_index': 5},
                          {'res1_index': 3, 'res2_index': 4}]
        ppv = c._ppv_score(input_contacts, RCmap, RCmap_seq)
        self.assertEqual(0.4, ppv)

    def test_refContacts(self):
        input = [[19.089, 51.837, 57.937], [16.279, 47.215, 58.769],
                 [7.599, 48.942, 52.633], [14.616, 50.231, 60.890],
                 [6.783, 46.760, 42.296]]
        input = [numpy.array(x) for x in input]

        length = 5
        dist_mat = numpy.zeros((length, length), numpy.float)
        dist_mat[0,0] = numpy.sqrt(numpy.sum((input[0] - input[0]) ** 2))
        dist_mat[0,1] = numpy.sqrt(numpy.sum((input[0] - input[1]) ** 2))
        dist_mat[0,2] = numpy.sqrt(numpy.sum((input[0] - input[2]) ** 2))
        dist_mat[0,3] = numpy.sqrt(numpy.sum((input[0] - input[3]) ** 2))
        dist_mat[0,4] = numpy.sqrt(numpy.sum((input[0] - input[4]) ** 2))
        dist_mat[1,0] = numpy.sqrt(numpy.sum((input[1] - input[0]) ** 2))
        dist_mat[1,1] = numpy.sqrt(numpy.sum((input[1] - input[1]) ** 2))
        dist_mat[1,2] = numpy.sqrt(numpy.sum((input[1] - input[2]) ** 2))
        dist_mat[1,3] = numpy.sqrt(numpy.sum((input[1] - input[3]) ** 2))
        dist_mat[1,4] = numpy.sqrt(numpy.sum((input[1] - input[4]) ** 2))
        dist_mat[2,0] = numpy.sqrt(numpy.sum((input[2] - input[0]) ** 2))
        dist_mat[2,1] = numpy.sqrt(numpy.sum((input[2] - input[1]) ** 2))
        dist_mat[2,2] = numpy.sqrt(numpy.sum((input[2] - input[2]) ** 2))
        dist_mat[2,3] = numpy.sqrt(numpy.sum((input[2] - input[3]) ** 2))
        dist_mat[2,4] = numpy.sqrt(numpy.sum((input[2] - input[4]) ** 2))
        dist_mat[3,0] = numpy.sqrt(numpy.sum((input[3] - input[0]) ** 2))
        dist_mat[3,1] = numpy.sqrt(numpy.sum((input[3] - input[1]) ** 2))
        dist_mat[3,2] = numpy.sqrt(numpy.sum((input[3] - input[2]) ** 2))
        dist_mat[3,3] = numpy.sqrt(numpy.sum((input[3] - input[3]) ** 2))
        dist_mat[3,4] = numpy.sqrt(numpy.sum((input[3] - input[4]) ** 2))
        dist_mat[4,0] = numpy.sqrt(numpy.sum((input[4] - input[0]) ** 2))
        dist_mat[4,1] = numpy.sqrt(numpy.sum((input[4] - input[1]) ** 2))
        dist_mat[4,2] = numpy.sqrt(numpy.sum((input[4] - input[2]) ** 2))
        dist_mat[4,3] = numpy.sqrt(numpy.sum((input[4] - input[3]) ** 2))
        dist_mat[4,4] = numpy.sqrt(numpy.sum((input[4] - input[4]) ** 2))

        ref_map = dist_mat < 8
        ref_contacts = numpy.where(dist_mat < 8)
        
        c = contacts_util.Contacter()
        contacts_out, contact_map = c._cb_contacts(input, input, length)
        
        numpy.testing.assert_equal(ref_contacts, contacts_out)
        numpy.testing.assert_equal(ref_map, contact_map)

    def test_truncate(self):
        input = [{'res1_index': 1, 'res2_index': 10},
                 {'res1_index': 1, 'res2_index': 6},
                 {'res1_index': 1, 'res2_index': 5},
                 {'res1_index': 1, 'res2_index': 2}]
        
        c = contacts_util.Contacter()
        ref_c = [{'res1_index': 1, 'res2_index': 10},
                 {'res1_index': 1, 'res2_index': 6},
                 {'res1_index': 1, 'res2_index': 5},
                 {'res1_index': 1, 'res2_index': 2}]
        out_c = c._truncateContactList(input, int(4*1.0))
        self.assertEqual(ref_c, out_c)
         
        ref_c = [{'res1_index': 1, 'res2_index': 10},
                 {'res1_index': 1, 'res2_index': 6}]
        out_c = c._truncateContactList(input, int(4*0.7))
        self.assertEqual(ref_c, out_c)
        
        ref_c = []
        out_c = c._truncateContactList(input, int(4*0.0))
        self.assertEqual(ref_c, out_c)
        
        ref_c = [{'res1_index': 1, 'res2_index': 10},
                 {'res1_index': 1, 'res2_index': 6},
                 {'res1_index': 1, 'res2_index': 5},
                 {'res1_index': 1, 'res2_index': 2}]
        out_c = c._truncateContactList(input, int(4*2.0))
        self.assertEqual(ref_c, out_c)

if __name__ == "__main__":
    unittest.main()    
