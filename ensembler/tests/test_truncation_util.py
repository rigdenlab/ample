"""Test functions for ensembler.truncation_util.py"""

import collections
import unittest

from ample.ensembler import truncation_util 

class Test(unittest.TestCase):
    
    def testResiduesFocussed(self):

        TheseusVariances = collections.namedtuple('TheseusVariances', 
                                                  ['idx', 'resName', 'resSeq', 
                                                   'variance', 'stdDev', 'rmsd', 'core'])
        l = 160
        var_by_res = [ TheseusVariances(idx=i, 
                                        resName='', 
                                        resSeq=i, 
                                        variance=float(i+1), 
                                        stdDev=None, rmsd=None, core=None) for i in range(l) ]

        ref_tlevels = [100, 93, 85, 78, 70, 63, 55, 48, 40, 33, 25, 23, 20, 18, 
                       15, 13, 10, 8, 5, 3]
        ref_tvariances = [160.0, 148.0, 136.0, 124.0, 112.0, 100.0, 88.0, 76.0, 
                          64.0, 52.0, 40.0, 36.0, 32.0, 28.0, 24.0, 20.0, 16.0, 
                          12.0, 8.0, 4.0]
        ref_tresidues = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 
                         [0, 1, 2, 3, 4, 5, 6, 7], 
                         [0, 1, 2, 3]]
        ref_tresidue_idxs = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 
                             [0, 1, 2, 3, 4, 5, 6, 7], 
                             [0, 1, 2, 3]]

        tlevels, tvariances, tresidues, tresidue_idxs = truncation_util.calculate_residues_focussed(var_by_res)
         
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tvariances, tvariances)
        for ref_tresidue, tresidue in zip(ref_tresidues, tresidues):
            self.assertEqual(ref_tresidue, tresidue)
        for ref_tresidue_idx, tresidue_idx in zip(ref_tresidue_idxs, tresidue_idxs):
            self.assertEqual(ref_tresidue_idx, tresidue_idx)
            
        return
    
    def testResiduesPercent(self):

        TheseusVariances = collections.namedtuple("TheseusVariances", 
                                                  ["idx", "resName", "resSeq", 
                                                   "variance", "stdDev", "rmsd", "core"])
        
        var_by_res = [TheseusVariances(idx=0, resName='GLN', resSeq=1, variance=55.757579, stdDev=7.4671, rmsd=18.603266, core=True), 
                      TheseusVariances(idx=1, resName='PRO', resSeq=2, variance=46.981224, stdDev=6.854285, rmsd=17.076522, core=True), 
                      TheseusVariances(idx=2, resName='ARG', resSeq=3, variance=47.734219, stdDev=6.908996, rmsd=17.212826, core=True), 
                      TheseusVariances(idx=3, resName='ARG', resSeq=4, variance=39.857312, stdDev=6.313265, rmsd=15.728643, core=True), 
                      TheseusVariances(idx=4, resName='LYS', resSeq=5, variance=35.477419, stdDev=5.956292, rmsd=14.839295, core=True), 
                      TheseusVariances(idx=5, resName='LEU', resSeq=6, variance=26.066709, stdDev=5.105557, rmsd=12.719802, core=True), 
                      TheseusVariances(idx=6, resName='CYS', resSeq=7, variance=24.114483, stdDev=4.91065, rmsd=12.234218, core=True), 
                      TheseusVariances(idx=7, resName='ILE', resSeq=8, variance=24.610979, stdDev=4.960945, rmsd=12.359523, core=True), 
                      TheseusVariances(idx=8, resName='LEU', resSeq=9, variance=21.187131, stdDev=4.602948, rmsd=11.467621, core=True), 
                      TheseusVariances(idx=9, resName='HIS', resSeq=10, variance=21.882362, stdDev=4.677859, rmsd=11.654251, core=True), 
                      TheseusVariances(idx=10, resName='ARG', resSeq=11, variance=21.62225, stdDev=4.649973, rmsd=11.584778, core=True), 
                      TheseusVariances(idx=11, resName='ASN', resSeq=12, variance=18.680585, stdDev=4.322104, rmsd=10.767937, core=True), 
                      TheseusVariances(idx=12, resName='PRO', resSeq=13, variance=16.568056, stdDev=4.070388, rmsd=10.140819, core=True), 
                      TheseusVariances(idx=13, resName='GLY', resSeq=14, variance=14.889562, stdDev=3.858699, rmsd=9.613426, core=True),
                      TheseusVariances(idx=14, resName='ARG', resSeq=15, variance=13.889743, stdDev=3.726895, rmsd=9.285052, core=True), 
                      TheseusVariances(idx=15, resName='CYS', resSeq=16, variance=8.722879, stdDev=2.953452, rmsd=7.358125, core=True), 
                      TheseusVariances(idx=16, resName='TYR', resSeq=17, variance=8.719477, stdDev=2.952876, rmsd=7.35669, core=True), 
                      TheseusVariances(idx=17, resName='ASP', resSeq=18, variance=4.648089, stdDev=2.155943, rmsd=5.371239, core=True), 
                      TheseusVariances(idx=18, resName='LYS', resSeq=19, variance=4.263944, stdDev=2.064932, rmsd=5.144498, core=True), 
                      TheseusVariances(idx=19, resName='ILE', resSeq=20, variance=2.338536, stdDev=1.529227, rmsd=3.809862, core=True)]
     
        # Test Case 1
        ref_tlevels = [100, 95, 90, 85, 80, 75, 70, 65, 60, 
                       55, 50, 45, 40, 35, 30, 25, 20, 15]
        ref_tvariances = [55.757579, 47.734219, 46.981224, 39.857312, 
                          35.477419, 26.066709, 24.610979, 24.114483, 
                          21.882362, 21.62225, 21.187131, 18.680585, 
                          16.568056, 14.889562, 13.889743, 8.722879, 
                          8.719477, 4.648089]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [9, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [13, 14, 15, 16, 17, 18, 19, 20], 
                         [14, 15, 16, 17, 18, 19, 20], 
                         [15, 16, 17, 18, 19, 20], 
                         [16, 17, 18, 19, 20], 
                         [17, 18, 19, 20], 
                         [18, 19, 20]]
        tlevels, tvariances, tresidues, _ = truncation_util.calculate_residues_percent(var_by_res, percent_interval=5)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        self.assertEqual(ref_tvariances, tvariances)
        
        ########################################################################       
        # Test Case 2
        ref_tlevels = [100, 85, 70, 55, 40, 25]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [13, 14, 15, 16, 17, 18, 19, 20], 
                         [16, 17, 18, 19, 20]]
        ref_tvariances = [55.757579, 39.857312, 24.610979, 21.62225, 
                          16.568056, 8.722879]
        tlevels, tvariances, tresidues, _ = truncation_util.calculate_residues_percent(var_by_res, percent_interval=15)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        self.assertEqual(ref_tvariances, tvariances)
        
        ########################################################################       
        # Test Case 3
        ref_tlevels = [100, 50]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 
                         [9, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
        ref_tvariances = [55.757579, 21.187131]
        tlevels, tvariances, tresidues, _ = truncation_util.calculate_residues_percent(var_by_res, percent_interval=50)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        self.assertEqual(ref_tvariances, tvariances)
        
        return
    
    def testResiduesThresh(self):
        '''Test we can calculate the original list of residues'''
        
        TheseusVariances = collections.namedtuple("TheseusVariances", 
                                                  ["idx", "resName", "resSeq", 
                                                   "variance", "stdDev", "rmsd", "core"])
        
        var_by_res = [TheseusVariances(idx=0, resName='GLN', resSeq=1, variance=55.757579, stdDev=7.4671, rmsd=18.603266, core=True), 
                      TheseusVariances(idx=1, resName='PRO', resSeq=2, variance=46.981224, stdDev=6.854285, rmsd=17.076522, core=True), 
                      TheseusVariances(idx=2, resName='ARG', resSeq=3, variance=47.734219, stdDev=6.908996, rmsd=17.212826, core=True), 
                      TheseusVariances(idx=3, resName='ARG', resSeq=4, variance=39.857312, stdDev=6.313265, rmsd=15.728643, core=True), 
                      TheseusVariances(idx=4, resName='LYS', resSeq=5, variance=35.477419, stdDev=5.956292, rmsd=14.839295, core=True), 
                      TheseusVariances(idx=5, resName='LEU', resSeq=6, variance=26.066709, stdDev=5.105557, rmsd=12.719802, core=True), 
                      TheseusVariances(idx=6, resName='CYS', resSeq=7, variance=24.114483, stdDev=4.91065, rmsd=12.234218, core=True), 
                      TheseusVariances(idx=7, resName='ILE', resSeq=8, variance=24.610979, stdDev=4.960945, rmsd=12.359523, core=True), 
                      TheseusVariances(idx=8, resName='LEU', resSeq=9, variance=21.187131, stdDev=4.602948, rmsd=11.467621, core=True), 
                      TheseusVariances(idx=9, resName='HIS', resSeq=10, variance=21.882362, stdDev=4.677859, rmsd=11.654251, core=True), 
                      TheseusVariances(idx=10, resName='ARG', resSeq=11, variance=21.62225, stdDev=4.649973, rmsd=11.584778, core=True), 
                      TheseusVariances(idx=11, resName='ASN', resSeq=12, variance=18.680585, stdDev=4.322104, rmsd=10.767937, core=True), 
                      TheseusVariances(idx=12, resName='PRO', resSeq=13, variance=16.568056, stdDev=4.070388, rmsd=10.140819, core=True), 
                      TheseusVariances(idx=13, resName='GLY', resSeq=14, variance=14.889562, stdDev=3.858699, rmsd=9.613426, core=True),
                      TheseusVariances(idx=14, resName='ARG', resSeq=15, variance=13.889743, stdDev=3.726895, rmsd=9.285052, core=True)]

        ########################################################################       
        # Test Case 1
        ref_tlevels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [7, 9, 10, 11, 12, 13, 14, 15], 
                         [9, 10, 11, 12, 13, 14, 15], 
                         [9, 11, 12, 13, 14, 15], 
                         [9, 12, 13, 14, 15], 
                         [12, 13, 14, 15], 
                         [13, 14, 15], 
                         [14, 15], 
                         [15]]
        tlevels, _, tresidues, _ = truncation_util.calculate_residues_thresh(var_by_res, percent_interval=10)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        
        ########################################################################       
        # Test Case 2
        ref_tlevels = [1, 2, 3, 4, 5]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [9, 11, 12, 13, 14, 15], 
                         [13, 14, 15]]
        
        tlevels, _, tresidues, _ = truncation_util.calculate_residues_thresh(var_by_res, percent_interval=20)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        
        ########################################################################       
        # Test Case 3
        ref_tlevels = [1, 2, 3]
        ref_tresidues = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
                         [9, 10, 11, 12, 13, 14, 15]]
        tlevels, _, tresidues, _ = truncation_util.calculate_residues_thresh(var_by_res, percent_interval=50)
        self.assertEqual(ref_tlevels, tlevels)
        self.assertEqual(ref_tresidues, tresidues)
        
        return
   
    def testGenerateThresholds(self):
        '''Test for Jaclyn's threshold method'''

        var_list = [31.186, 24.104, 16.448, 40.365, 52.0, 11.549999999999999, 
                    18.252000000000002, 9.324000000000002, 19.551, 46.53, 
                    21.609, 1.0090000000000001, 14.504, 13.39, 7.259, 33.231, 
                    29.54, 20.64, 25.325000000000003, 34.947, 7.147, 50.85, 
                    44.290000000000006, 8.312, 58.254000000000005, 
                    53.907000000000004, 50.1, 2.028, 31.65, 50.666000000000004]
        
        # Test Case 1
        ref_thresholds = [1.0090000000000001, 2.028, 7.147, 7.259, 8.312, 
                          9.324000000000002, 11.549999999999999, 13.39, 
                          14.504, 16.448, 18.252000000000002, 19.551, 20.64, 
                          21.609, 24.104, 25.325000000000003, 29.54, 31.186, 
                          31.65, 33.231, 34.947, 40.365, 44.290000000000006, 
                          46.53, 50.1, 50.666000000000004, 50.85, 52.0,
                          53.907000000000004, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(5))
        thresholds = truncation_util._generate_thresholds(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)
        
        # Test Case 2
        ref_thresholds = [9.324000000000002, 19.551, 31.186, 46.53, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(20))
        thresholds = truncation_util._generate_thresholds(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)
        
        # Test Case 3
        ref_thresholds = [24.104, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(50))
        thresholds = truncation_util._generate_thresholds(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)
        
        return
    
    def testGenerateThresholds2(self):
        '''Test for Jens' threshold method'''

        var_list = [31.186, 24.104, 16.448, 40.365, 52.0, 11.549999999999999, 
                    18.252000000000002, 9.324000000000002, 19.551, 46.53, 
                    21.609, 1.0090000000000001, 14.504, 13.39, 7.259, 33.231, 
                    29.54, 20.64, 25.325000000000003, 34.947, 7.147, 50.85, 
                    44.290000000000006, 8.312, 58.254000000000005, 
                    53.907000000000004, 50.1, 2.028, 31.65, 50.666000000000004]
        
        # Test Case 1
        ref_thresholds = [1.0090000000000001, 2.028, 7.147, 7.259, 8.312, 
                          9.324000000000002, 11.549999999999999, 13.39, 
                          14.504, 16.448, 18.252000000000002, 19.551, 20.64, 
                          21.609, 24.104, 25.325000000000003, 29.54, 31.186, 
                          31.65, 33.231, 34.947, 40.365, 44.290000000000006, 
                          46.53, 50.1, 50.666000000000004, 50.85, 52.0, 
                          53.907000000000004, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(5))
        thresholds = truncation_util._generate_thresholds2(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)

        # Test Case 2
        ref_thresholds = [9.324000000000002, 19.551, 31.186, 46.53, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(20))
        thresholds = truncation_util._generate_thresholds2(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)

        # Test Case 3
        ref_thresholds = [24.104, 58.254000000000005]
        chunk_size = int((float(len(var_list)) / 100) * float(50))
        thresholds = truncation_util._generate_thresholds2(var_list, chunk_size)
        self.assertEqual(ref_thresholds, thresholds)
        
        return  
    
    def testPruneResidues(self):
        '''Test we can reproduce the original thresholds'''

        residues = []
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, None)
          
        residues = [1]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, [1])
         
        residues = [1, 2]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [])
        self.assertEqual(pruned, [1, 2])
         
        residues = [1, 2, 4]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 4])
        self.assertEqual(pruned, None)
         
        # Big enough gap 
        residues = [1, 2, 3, 4, 8, 13, 14, 15]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 3, 4, 13, 14, 15])
        self.assertEqual(pruned, [8])
           
        residues = [1, 2, 3, 4, 8, 9, 13, 14, 15]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [1, 2, 3, 4, 13, 14, 15])
        self.assertEqual(pruned, [8, 9])
           
        residues = [1, 2, 3, 4, 8, 9, 10, 13, 14, 15]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, residues)
        self.assertEqual(pruned, None)
           
        # end gap not big enough
        residues = [1, 2, 3, 4, 8, 9, 11, 12, 13]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, residues)
        self.assertEqual(pruned, None)
           
        # Lone residue at start
        residues = [1, 11, 12, 13]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [11, 12, 13])
        self.assertEqual(pruned, [1])
        
        # Lone residue at end
        residues = [11, 12, 13, 19]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres, [11, 12, 13])
        self.assertEqual(pruned, [19])
          
        # Mixed
        residues = [1, 3, 4, 7, 10, 11, 13, 15, 16, 19]
        pres, pruned = truncation_util.prune_residues(residues, chunk_size=1, allowed_gap=2)
        self.assertEqual(pres, [1, 3, 4, 10, 11, 13, 15, 16])
        self.assertEqual(pruned, [7, 19])

        return
    
    def testSplitSequence(self):
        
        ########################################################################       
        # Test Case 1
        ref_idxs = [99, 94, 89, 84, 79, 74, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 9, 4]
        idxs = truncation_util._split_sequence(100, 5, 3)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 2        
        ref_idxs = [99, 79, 59, 39, 19]
        idxs = truncation_util._split_sequence(100, 20, 3)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 3        
        ref_idxs = [99, 49]
        idxs = truncation_util._split_sequence(100, 50, 3)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 4        
        ref_idxs = [99, 94, 89, 84, 79, 74, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24]
        idxs = truncation_util._split_sequence(100, 5, 21)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 5        
        ref_idxs = [99, 79, 59, 39]
        idxs = truncation_util._split_sequence(100, 20, 21)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 6        
        ref_idxs = [99, 49]
        idxs = truncation_util._split_sequence(100, 50, 21)
        self.assertEqual(ref_idxs, idxs)
        
        ########################################################################       
        # Test Case 7
        ref_idxs = [99]
        idxs = truncation_util._split_sequence(100, 50, 99)
        self.assertEqual(ref_idxs, idxs)
                
        return
    
if __name__ == "__main__":
    unittest.main()
