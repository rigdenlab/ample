"""Test facility for ample/util/pdb_edit.py"""

__author__ = "Jens Thomas, Adam Simpkin & Felix Simkovic"
__date__ = "22 Apr 2017"

import glob
import iotbx.pdb
import logging
import os
import tempfile
import unittest

from ample import constants
from ample.util import pdb_edit


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[:-2])
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(os.sep.join(paths[:-3]), 'testfiles')

    def test_backbone_1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp_4DZN.pdb"
        reference_data = [
            "ATOM      4  N   GLY A   1      24.076  12.643  -9.179  1.00 23.32           N",
            "ATOM      5  CA  GLY A   1      22.806  12.124  -9.698  1.00 18.13           C",
            "ATOM      6  C   GLY A   1      22.170  11.067  -8.799  1.00 15.67           C",
            "ATOM      7  O   GLY A   1      22.404  11.024  -7.580  1.00 16.52           O",
            "ATOM      8  N   GLU A   2      21.377  10.190  -9.397  1.00 13.90           N",
            "ATOM      9  CA  GLU A   2      20.675   9.156  -8.637  1.00 11.97           C",
            "ATOM     10  C   GLU A   2      21.614   8.106  -7.996  1.00 13.90           C",
            "ATOM     11  O   GLU A   2      21.337   7.619  -6.898  1.00 12.88           O",
            "ATOM     12  CB  GLU A   2      19.625   8.485  -9.531  1.00 13.40           C",
            "ATOM     17  N   ILE A   3      22.722   7.760  -8.656  1.00 10.53           N",
        ]
        pdb_edit.backbone(pdbin, pdbout)
        data = [line[0:78] for line in open(pdbout) if line.startswith("ATOM")]
        self.assertEqual(data[:10], reference_data)
        os.unlink(pdbout)

    def test_calpha_only_1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp_4DZN.pdb"
        reference_data = [
            "ATOM      5  CA  GLY A   1      22.806  12.124  -9.698  1.00 18.13           C",
            "ATOM      9  CA  GLU A   2      20.675   9.156  -8.637  1.00 11.97           C",
            "ATOM     18  CA  ILE A   3      23.668   6.824  -8.067  1.00 14.88           C",
            "ATOM     28  CA  ALA A   4      25.130   9.427  -5.671  1.00 16.59           C",
            "ATOM     33  CA  ALA A   5      21.738   9.690  -3.986  1.00 12.84           C",
            "ATOM     38  CA  LEU A   6      21.616   5.898  -3.440  1.00 10.59           C",
            "ATOM     46  CA  LYS A   7      25.138   6.027  -1.978  1.00 11.22           C",
            "ATOM     55  CA  GLN A   8      23.934   8.766   0.443  1.00 12.72           C",
            "ATOM     64  CA  GLU A   9      21.017   6.487   1.401  1.00 12.29           C",
            "ATOM     73  CA  ILE A  10      23.439   3.652   2.166  1.00 10.77           C",
        ]
        pdb_edit.calpha_only(pdbin, pdbout)
        data = [line[0:78] for line in open(pdbout) if line.startswith("ATOM")]
        self.assertEqual(data[:10], reference_data)
        os.unlink(pdbout)

    def test_extract_chain_1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp_4DZN.pdb"
        reference_data = [
            "ATOM      1  N   GLY B   1      24.076  12.643  -9.179  1.00 23.32           N",
            "ATOM      2  CA  GLY B   1      22.806  12.124  -9.698  1.00 18.13           C",
            "ATOM      3  C   GLY B   1      22.170  11.067  -8.799  1.00 15.67           C",
            "ATOM      4  O   GLY B   1      22.404  11.024  -7.580  1.00 16.52           O",
            "ATOM      5  N   GLU B   2      21.377  10.190  -9.397  1.00 13.90           N",
            "ATOM      6  CA  GLU B   2      20.675   9.156  -8.637  1.00 11.97           C",
            "ATOM      7  C   GLU B   2      21.614   8.106  -7.996  1.00 13.90           C"
        ]
        pdb_edit.extract_chain(pdbin, pdbout, 'A', new_chain_id='B', c_alpha=False, renumber=True)
        data = [line[0:78] for line in open(pdbout) if line.startswith("ATOM")]
        self.assertEqual(data[:7], reference_data)
        os.unlink(pdbout)

    def test_extract_chain_2(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp_4DZN.pdb"
        reference_data = [
            "ATOM      1  CA  GLY A   1      22.806  12.124  -9.698  1.00 18.13           C",
            "ATOM      2  CA  GLU A   2      20.675   9.156  -8.637  1.00 11.97           C",
            "ATOM      3  CA  ILE A   3      23.668   6.824  -8.067  1.00 14.88           C",
            "ATOM      4  CA  ALA A   4      25.130   9.427  -5.671  1.00 16.59           C",
            "ATOM      5  CA  ALA A   5      21.738   9.690  -3.986  1.00 12.84           C",
            "ATOM      6  CA  LEU A   6      21.616   5.898  -3.440  1.00 10.59           C",
            "ATOM      7  CA  LYS A   7      25.138   6.027  -1.978  1.00 11.22           C"
        ]
        pdb_edit.extract_chain(pdbin, pdbout, 'A', new_chain_id=None, c_alpha=True, renumber=True)
        data = [line[0:78] for line in open(pdbout) if line.startswith("ATOM")]
        self.assertEqual(data[:7], reference_data)
        os.unlink(pdbout)

    def test_extract_model_1(self):
        pdbin = os.path.join(self.ample_share, "examples", "import-data", "input",
                             "ensembles", "ensemble_1.pdb")
        pdbout = "tmp_ensemble_1.pdb"
        model_id = 1

        reference_data = [
            "ATOM      1  N   GLY A  15       9.003  -8.476  -2.672  1.00  0.00           N",
            "ATOM      2  CA  GLY A  15       7.637  -8.931  -2.443  1.00  0.00           C",
            "ATOM      3  C   GLY A  15       6.920  -8.037  -1.440  1.00  0.00           C",
            "ATOM      4  O   GLY A  15       5.871  -7.469  -1.740  1.00  0.00           O",
            "ATOM      5  H   GLY A  15       9.692  -9.043  -2.553  1.00  0.00           H",
            "ATOM      6  N   PRO A  16       7.494  -7.916  -0.247  1.00  0.00           N",
            "ATOM      7  CA  PRO A  16       6.889  -7.124   0.817  1.00  0.00           C",
            "ATOM      8  C   PRO A  16       6.560  -5.717   0.337  1.00  0.00           C",
            "ATOM      9  O   PRO A  16       5.559  -5.128   0.749  1.00  0.00           O",
            "ATOM     10  CB  PRO A  16       7.953  -7.106   1.917  1.00  0.00           C"
        ]
        pdb_edit.extract_model(pdbin, pdbout, model_id)
        data = [line[0:78] for line in open(pdbout) if line.startswith("ATOM")]
        self.assertEqual(data[:10], reference_data)
        os.unlink(pdbout)

    def test_extract_resSeq_1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        data = pdb_edit.extract_resSeq(pdbin)
        ref_data = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                    30, 31, 32]
        self.assertEqual(data, ref_data)

    def test_extract_resSeq_2(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        data = pdb_edit.extract_resSeq(pdbin, chain_id='B')
        ref_data = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                    30, 31, 32]
        self.assertEqual(data, ref_data)

    def test_keep_residues_1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp_4DZN.pdb"
        residue_range = [5, 10]
        chain_id = "A"
        reference_data = [
            "ATOM     27  N   ALA A   5      24.623   8.709  -6.856  1.00 14.20           N",
            "ATOM     28  CA  ALA A   5      25.130   9.427  -5.671  1.00 16.59           C",
            "ATOM     29  C   ALA A   5      24.079   9.482  -4.565  1.00 13.45           C",
            "ATOM     30  O   ALA A   5      24.399   9.349  -3.366  1.00 14.72           O",
            "ATOM     31  CB  ALA A   5      25.599  10.846  -6.036  1.00 18.01           C",
            "ATOM     32  N   ALA A   6      22.825   9.671  -4.958  1.00 11.84           N",
            "ATOM     33  CA  ALA A   6      21.738   9.690  -3.986  1.00 12.84           C",
            "ATOM     34  C   ALA A   6      21.566   8.342  -3.261  1.00 15.19           C",
            "ATOM     35  O   ALA A   6      21.356   8.313  -2.061  1.00 13.05           O",
            "ATOM     36  CB  ALA A   6      20.430  10.143  -4.644  1.00 14.15           C",
            "ATOM     37  N   LEU A   7      21.667   7.244  -4.001  1.00 12.44           N",
            "ATOM     38  CA  LEU A   7      21.616   5.898  -3.440  1.00 10.59           C",
            "ATOM     39  C   LEU A   7      22.784   5.647  -2.489  1.00 11.89           C",
            "ATOM     40  O   LEU A   7      22.601   5.077  -1.428  1.00 12.64           O",
            "ATOM     41  CB  LEU A   7      21.570   4.834  -4.538  1.00 13.99           C",
            "ATOM     42  CG  LEU A   7      20.240   4.827  -5.307  1.00 14.68           C",
            "ATOM     43  CD1 LEU A   7      20.342   4.020  -6.589  1.00 20.23           C",
            "ATOM     44  CD2 LEU A   7      19.130   4.307  -4.425  1.00 21.58           C",
            "ATOM     45  N   LYS A   8      23.970   6.104  -2.865  1.00  8.95           N",
            "ATOM     46  CA  LYS A   8      25.138   6.027  -1.978  1.00 11.22           C",
            "ATOM     47  C   LYS A   8      24.909   6.801  -0.665  1.00 13.03           C",
            "ATOM     48  O   LYS A   8      25.275   6.317   0.399  1.00 13.72           O",
            "ATOM     49  CB  LYS A   8      26.376   6.546  -2.680  1.00 14.45           C",
            "ATOM     50  CG  LYS A   8      26.897   5.631  -3.763  1.00 16.58           C",
            "ATOM     51  CD  LYS A   8      28.122   6.277  -4.424  1.00 21.80           C",
            "ATOM     52  CE  LYS A   8      28.803   5.351  -5.369  0.50 24.13           C",
            "ATOM     53  NZ  LYS A   8      29.986   6.016  -5.957  0.50 18.04           N",
            "ATOM     54  N   GLN A   9      24.313   7.988  -0.754  1.00 12.77           N",
            "ATOM     55  CA  GLN A   9      23.934   8.766   0.443  1.00 12.72           C",
            "ATOM     56  C   GLN A   9      22.917   8.009   1.310  1.00 12.01           C",
            "ATOM     57  O   GLN A   9      22.998   8.034   2.553  1.00 14.90           O",
            "ATOM     58  CB  GLN A   9      23.402  10.160   0.080  1.00 22.08           C",
            "ATOM     59  CG  GLN A   9      23.223  11.191   0.978  0.00 34.62           C",
            "ATOM     60  CD  GLN A   9      22.739  12.519   0.396  0.00 45.70           C",
            "ATOM     61  OE1 GLN A   9      22.841  12.760  -0.810  0.00 49.00           O",
            "ATOM     62  NE2 GLN A   9      22.210  13.385   1.258  0.00 44.32           N",
            "ATOM     63  N   GLU A  10      21.979   7.310   0.671  1.00 11.84           N",
            "ATOM     64  CA  GLU A  10      21.017   6.487   1.401  1.00 12.29           C",
            "ATOM     65  C   GLU A  10      21.717   5.383   2.184  1.00 10.63           C",
            "ATOM     66  O   GLU A  10      21.390   5.142   3.344  1.00 13.04           O",
            "ATOM     67  CB  GLU A  10      19.974   5.843   0.475  1.00 17.18           C",
            "ATOM     68  CG  GLU A  10      19.043   6.797  -0.218  1.00 21.50           C",
            "ATOM     69  CD  GLU A  10      17.852   7.252   0.598  0.50 26.52           C",
            "ATOM     70  OE1 GLU A  10      17.710   6.852   1.763  0.50 21.11           O",
            "ATOM     71  OE2 GLU A  10      17.047   8.032   0.047  0.50 28.78           O"
        ]
        pdb_edit.keep_residues(pdbin, pdbout, residue_range, chain_id)
        lines_to_read = range(6, 51)
        data = []
        with open(pdbout) as f:
            for i, line in enumerate(f):
                if i in lines_to_read:
                    data.append(line[0:78])
                else:
                    continue
        self.assertEqual(data, reference_data)
        os.unlink(pdbout)

    def test_merge(self):
        pdb1 = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdb2 = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = "tmp_merge.pdb"
        pdb_edit.merge(pdb1, pdb2, pdbout)
        reference_data_1 = [
            "HETATM    1  C   ACE A   0      25.199  11.913  -9.250  1.00 27.72           C",
            "ANISOU    1  C   ACE A   0     2933   4198   3402     29   -251    297       C",
            "HETATM    2  O   ACE A   0      25.201  10.666  -9.372  1.00 23.97           O",
            "ANISOU    2  O   ACE A   0     2587   4332   2189   -515   -230    104       O",
            "HETATM    3  CH3 ACE A   0      26.454  12.702  -9.001  1.00 32.42           C",
            "ANISOU    3  CH3 ACE A   0     3564   4758   3995   -239  -1016     28       C",
            "ATOM      4  N   GLY A   1      24.076  12.643  -9.179  1.00 23.32           N",
            "ANISOU    4  N   GLY A   1     2998   3289   2573   -157   -425    199       N",
            "ATOM      5  CA  GLY A   1      22.806  12.124  -9.698  1.00 18.13           C",
            "ANISOU    5  CA  GLY A   1     2466   2494   1928    -66    -40    660       C",
            "ATOM      6  C   GLY A   1      22.170  11.067  -8.799  1.00 15.67           C",
            "ANISOU    6  C   GLY A   1     2501   1562   1889   -506   -480    259       C",
            "ATOM      7  O   GLY A   1      22.404  11.024  -7.580  1.00 16.52           O",
            "ANISOU    7  O   GLY A   1     2581   2124   1571   -470   -178    -48       O",
            "ATOM      8  N   GLU A   2      21.377  10.190  -9.397  1.00 13.90           N",
        ]

        reference_data_2 = [
            "ANISOU  576  CB  ILE C   3     2180   1836   1869   -499   -391    255       C",
            "ATOM    577  CG1 ILE C   3      19.687   0.791 -11.270  1.00 16.36           C",
            "ANISOU  577  CG1 ILE C   3     2342   1902   1970    708   -222    402       C",
            "ATOM    578  CG2 ILE C   3      17.926   2.092 -10.024  1.00 16.72           C",
            "ANISOU  578  CG2 ILE C   3     2270   2100   1980    393    314    -53       C",
            "ATOM    579  CD1 ILE C   3      19.861   1.505 -12.619  1.00 19.26           C",
            "ANISOU  579  CD1 ILE C   3     2386   2732   2196     94   -787   1035       C",
            "ATOM    580  N   ALA C   4      15.647  -0.630  -9.580  1.00 13.96           N",
            "ANISOU  580  N   ALA C   4     1637   2255   1413   -554    -86    436       N",
            "ATOM    581  CA  ALA C   4      14.367  -0.707  -8.884  1.00 12.18           C",
            "ANISOU  581  CA  ALA C   4     1278   1976   1372   -216   -201   -146       C",
            "ATOM    582  C   ALA C   4      14.424  -1.670  -7.682  1.00 14.16           C",
            "ANISOU  582  C   ALA C   4     1516   1953   1911   -420   -221    -33       C",
            "ATOM    583  O   ALA C   4      13.931  -1.344  -6.603  1.00 15.41           O",
            "ANISOU  583  O   ALA C   4     1669   2382   1804   -498   -191   -137       O"
        ]

        lines_to_read_1 = range(351, 366)
        lines_to_read_2 = range(1500, 1515)
        data_1 = []
        data_2 = []
        with open(pdbout) as f:
            for i, line in enumerate(f):
                if i in lines_to_read_1:
                    data_1.append(line[0:78])
                elif i in lines_to_read_2:
                    data_2.append(line[0:78])
                else:
                    continue
        os.unlink(pdbout)
        self.assertEqual(data_1, reference_data_1)
        self.assertEqual(data_2, reference_data_2)
        return

    def test_most_prob(self):
        # Keep the most probably conformer
        pdbin = os.path.join(self.testfiles_dir, "2UUI.pdb")
        pdb_input = iotbx.pdb.pdb_input(pdbin)
        hierarchy = pdb_input.construct_hierarchy()

        pdb_edit.most_prob(hierarchy)

        # Residues which are conformers
        residues_to_check = [62, 82, 84]
        data = []

        # extract atom names for residues which were conformers
        for model in hierarchy.models():
            for chain in model.chains():
                for residue_group in chain.residue_groups():
                    if residue_group.resseq_as_int() in residues_to_check:
                        for atom_group in residue_group.atom_groups():
                            for atom in atom_group.atoms():
                                data.append(atom.name.strip())

        # Only one set of atoms returned if function works correctly
        reference_data = ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2',
                          'N', 'C', 'O', 'CA', 'CB', 'SG',
                          'N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2']

        self.assertEqual(data, reference_data)
        return

    def test_num_atoms_and_residues(self):
        # Extract the number of residues and atoms in a pdb
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        ref_natoms = 1711
        ref_nresidues = 93

        natoms, nresidues = pdb_edit.num_atoms_and_residues(pdbin)

        self.assertEqual(natoms, ref_natoms)
        self.assertEqual(nresidues, ref_nresidues)

        # Second test with first flag
        ref_natoms = 252
        ref_nresidues = 33

        natoms, nresidues = pdb_edit.num_atoms_and_residues(pdbin, first=True)

        self.assertEqual(natoms, ref_natoms)
        self.assertEqual(nresidues, ref_nresidues)

        return

    def test_rename_chains(self):
        # Test the function to rename chains
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp.pdb"
        reference_data = ["ATOM      4  N   GLY Z   1      24.076  12.643  -9.179  1.00 23.32           N",
                          "ATOM      5  CA  GLY Z   1      22.806  12.124  -9.698  1.00 18.13           C",
                          "ATOM      6  C   GLY Z   1      22.170  11.067  -8.799  1.00 15.67           C",
                          "ATOM      7  O   GLY Z   1      22.404  11.024  -7.580  1.00 16.52           O",
                          "ATOM      8  N   GLU Z   2      21.377  10.190  -9.397  1.00 13.90           N",
                          "ATOM      9  CA  GLU Z   2      20.675   9.156  -8.637  1.00 11.97           C",
                          "ATOM     10  C   GLU Z   2      21.614   8.106  -7.996  1.00 13.90           C",
                          "ATOM     11  O   GLU Z   2      21.337   7.619  -6.898  1.00 12.88           O",
                          "ATOM     12  CB  GLU Z   2      19.625   8.485  -9.531  1.00 13.40           C",
                          "ATOM     13  CG  GLU Z   2      18.637   7.595  -8.790  1.00 10.75           C"
                          ]
        pdb_edit.rename_chains(pdbin, pdbout, fromChain=['A', 'B'], toChain=['Z', 'Y'])

        count = 0
        data = []
        with open(pdbout) as f:
            for line in f:
                if count < 10:
                    if line.startswith("ATOM"):
                        data.append(line[0:78])
                        count += 1
                    else:
                        continue
        os.unlink(pdbout)
        self.assertEqual(data, reference_data)
        return

    def test_renumber(self):
        # Test the function to renumber residues
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "tmp.pdb"

        reference_data = ["ATOM      4  N   GLY A   6      24.076  12.643  -9.179  1.00 23.32           N",
                          "ATOM      5  CA  GLY A   6      22.806  12.124  -9.698  1.00 18.13           C",
                          "ATOM      6  C   GLY A   6      22.170  11.067  -8.799  1.00 15.67           C",
                          "ATOM      7  O   GLY A   6      22.404  11.024  -7.580  1.00 16.52           O",
                          "ATOM      8  N   GLU A   7      21.377  10.190  -9.397  1.00 13.90           N",
                          "ATOM      9  CA  GLU A   7      20.675   9.156  -8.637  1.00 11.97           C",
                          "ATOM     10  C   GLU A   7      21.614   8.106  -7.996  1.00 13.90           C"
                          ]
        pdb_edit.renumber_residues(pdbin, pdbout, start=5)

        count = 0
        data = []
        with open(pdbout) as f:
            for line in f:
                if count < 7:
                    if line.startswith("ATOM"):
                        data.append(line[0:78])
                        count += 1
                    else:
                        continue
        os.unlink(pdbout)
        self.assertEqual(data, reference_data)
        return

    def test_standardise(self):
        # Test standardisation of a pdb file

        #######################################################################
        # Test case 1 - Testing standardisation
        #######################################################################

        pdbin = tempfile.NamedTemporaryFile("w", suffix='.pdb', delete=False)

        pdbin.write("""HETATM    1  O   HOH A   0      25.199  11.913  -9.250  1.00 27.72           O
ANISOU    1  O   HOH A   0     2933   4198   3402     29   -251    297       O
HETATM    2  O   HOH A   0      25.201  10.666  -9.372  1.00 23.97           O
ANISOU    2  O   HOH A   0     2587   4332   2189   -515   -230    104       O
HETATM    3  O   HOH A   0      26.454  12.702  -9.001  1.00 32.42           O
ANISOU    3  O   HOH A   0     3564   4758   3995   -239  -1016     28       O
ATOM      4  N   URA A   1      24.076  12.643  -9.179  1.00 23.32           N
ANISOU    4  N   URA A   1     2998   3289   2573   -157   -425    199       N
ATOM      5  CA  URA A   1      22.806  12.124  -9.698  1.00 18.13           C
ANISOU    5  CA  URA A   1     2466   2494   1928    -66    -40    660       C
ATOM      6  C   URA A   1      22.170  11.067  -8.799  1.00 15.67           C
ANISOU    6  C   URA A   1     2501   1562   1889   -506   -480    259       C
ATOM      7  O   URA A   1      22.404  11.024  -7.580  1.00 16.52           O
ANISOU    7  O   URA A   1     2581   2124   1571   -470   -178    -48       O
ATOM      8  N   GLU A   2      21.377  10.190  -9.397  1.00 13.90           N
ANISOU    8  N   GLU A   2     2146   2122   1012   -480    -81     50       N
ATOM      9  CA  GLU A   2      20.675   9.156  -8.637  1.00 11.97           C
ANISOU    9  CA  GLU A   2     1678   1707   1163   -350    115     31       C
ATOM     10  C   GLU A   2      21.614   8.106  -7.996  1.00 13.90           C
ANISOU   10  C   GLU A   2     1930   1893   1457   -280    128     26       C
ATOM     11  O   GLU A   2      21.337   7.619  -6.898  1.00 12.88           O
ANISOU   11  O   GLU A   2     1404   2192   1297   -378   -136    183       O
ATOM     12  CB  GLU A   2      19.625   8.485  -9.531  1.00 13.40           C
ANISOU   12  CB  GLU A   2     1703   2108   1279   -610    157     41       C
ATOM     13  CG  GLU A   2      18.637   7.595  -8.790  1.00 10.75           C
ANISOU   13  CG  GLU A   2      859   1761   1463   -466   -737     10       C
ATOM     14  CD  GLU A   2      17.652   8.361  -7.951  1.00 15.12           C
ANISOU   14  CD  GLU A   2     2189   2375   1180    249   -108    115       C
ATOM     15  OE1 GLU A   2      17.724   9.603  -7.887  1.00 21.69           O
ANISOU   15  OE1 GLU A   2     2427   2488   3323   -101   -364    580       O
ATOM     16  OE2 GLU A   2      16.786   7.706  -7.365  1.00 25.52           O
ANISOU   16  OE2 GLU A   2     2703   4739   2252   -883    385   -267       O
ATOM     17  N   ILE A   3      22.722   7.760  -8.656  1.00 10.53           N
ANISOU   17  N   ILE A   3     1371   1474   1154   -336   -217    311       N
ATOM     18  CA  ILE A   3      23.668   6.824  -8.067  1.00 14.88           C
ANISOU   18  CA  ILE A   3     1766   2158   1730    207   -104    484       C
ATOM     19  C   ILE A   3      24.256   7.428  -6.800  1.00 12.88           C
ANISOU   19  C   ILE A   3     1122   2110   1661    314    171    525       C
ATOM     20  O   ILE A   3      24.339   6.741  -5.780  1.00 13.76           O
ANISOU   20  O   ILE A   3     2174   1880   1171    341   -223    635       O
ATOM     21  CB  ILE A   3      24.783   6.398  -9.051  1.00 15.77           C
ANISOU   21  CB  ILE A   3     2561   1806   1623    299    -22    275       C
ATOM     22  CG1AILE A   3      24.205   5.544 -10.195  0.50 24.17           C
ANISOU   22  CG1AILE A   3     4017   3103   2063    128   -307   -396       C
ATOM     23  CG1BILE A   3      24.158   5.509 -10.147  0.50 24.04           C
ANISOU   23  CG1BILE A   3     3933   3035   2166     67   -306   -400       C
ATOM     24  CG2 ILE A   3      25.906   5.657  -8.309  1.00 22.84           C
ANISOU   24  CG2 ILE A   3     3223   3912   1542   1636    151   -176       C
ATOM     25  CD1AILE A   3      23.657   4.220  -9.758  0.50 25.30           C
ANISOU   25  CD1AILE A   3     4137   3472   2001   -326    846   -761       C
ATOM     26  CD1BILE A   3      25.134   4.923 -11.133  0.50 25.34           C
ANISOU   26  CD1BILE A   3     4035   3461   2132   1269   -829   -356       C
""")
        pdbin.close()

        pdbout = "tmp_4DZN.pdb"
        reference_data = ["ATOM      1  N   GLU A   2      21.377  10.190  -9.397  1.00 13.90           N",
                          "ATOM      2  CA  GLU A   2      20.675   9.156  -8.637  1.00 11.97           C",
                          "ATOM      3  C   GLU A   2      21.614   8.106  -7.996  1.00 13.90           C",
                          "ATOM      4  O   GLU A   2      21.337   7.619  -6.898  1.00 12.88           O",
                          "ATOM      5  CB  GLU A   2      19.625   8.485  -9.531  1.00 13.40           C",
                          "ATOM      6  CG  GLU A   2      18.637   7.595  -8.790  1.00 10.75           C",
                          "ATOM      7  CD  GLU A   2      17.652   8.361  -7.951  1.00 15.12           C",
                          "ATOM      8  OE1 GLU A   2      17.724   9.603  -7.887  1.00 21.69           O",
                          "ATOM      9  OE2 GLU A   2      16.786   7.706  -7.365  1.00 25.52           O",
                          "ATOM     10  N   ILE A   3      22.722   7.760  -8.656  1.00 10.53           N"]

        pdb_edit.standardise(pdbin.name, pdbout, chain='A')

        count = 0
        data = []
        with open(pdbout) as f:
            for line in f:
                if count < 10:
                    # Make sure the water and anisou atoms are deleted
                    if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
                        data.append(line[0:78])
                        count += 1
                    else:
                        continue
        os.unlink(pdbout)
        self.assertEqual(data, reference_data)

    def test_std_Residues(self):

        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "std.pdb"

        pdb_edit.std_residues(pdbin, pdbout)

        # Check it's valid
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)

        # Get list of all the residue names in chain 1
        resnames = [g.unique_resnames()[0] for g in pdb_obj.hierarchy.models()[0].chains()[0].residue_groups()]
        ref = ['ACE', 'GLY', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'GLN', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU',
               'LYS', 'LYS', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'PHE', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU',
               'LYS', 'GLN', 'GLY', 'TYR', 'TYR']
        self.assertEqual(resnames, ref)

        os.unlink(pdbout)

        return

    def testGetInfo1(self):
        """"""

        pdbfile = os.path.join(self.testfiles_dir, "1GU8.pdb")

        info = pdb_edit.get_info(pdbfile)

        self.assertEqual(info.pdbCode, "1GU8")
        self.assertEqual(len(info.models), 2)

        m1 = info.models[0]
        self.assertEqual(m1.chains[0], 'A')
        self.assertEqual(m1.resSeqs[0],
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                          27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
                          73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                          96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
                          115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                          134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
                          153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
                          172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
                          191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
                          210, 211, 212, 213, 214, 215, 216, 217, 218, 219])
        self.assertEqual(m1.sequences[0],
                         'VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLD'
                         'SREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVAL'
                         'LTPTVDVALIVYLDLVTKVGFGFIALDAAATL')

        self.assertEqual(m1.caMask[0], [False] * 218)
        self.assertEqual(m1.bbMask[0],
                         [False, True, False, False, False, False, False, False, False, True, False, False, True, False,
                          False, False, True, False, False, False, False, False, False, False, True, False, False,
                          False, True, False, True, False, False, False, False, False, False, False, False, False, True,
                          False, False, True, False, False, False, False, False, False, False, False, False, False,
                          False, True, False, True, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, True, False, False, False, True, False, False,
                          False, False, False, False, True, False, False, False, False, False, False, False, False,
                          False, False, False, False, True, False, False, True, False, False, False, False, True, False,
                          False, False, False, False, False, False, True, False, True, False, False, False, False,
                          False, True, False, False, False, False, False, False, True, False, False, False, False,
                          False, False, False, False, False, False, False, True, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, True, False, False, True, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, True, False, True, False, False, False,
                          False, False, False, False, False, False, True])

        self.assertEqual(info.numAtoms(modelIdx=0), 1621)
        self.assertEqual(info.numCalpha(modelIdx=0), 218)

        m2 = info.models[1]
        self.assertEqual(m2.chains[0], 'A')
        self.assertEqual(m2.resSeqs[0],
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                          27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
                          73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                          96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
                          115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                          134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
                          153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
                          172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
                          191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
                          210, 211, 212, 213, 214, 215, 216, 217, 218, 219])
        self.assertEqual(m2.sequences[0],
                         'VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLD'
                         'SREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVAL'
                         'LTPTVDVALIVYLDLVTKVGFGFIALDAAATL')

        self.assertEqual(info.numAtoms(modelIdx=1), 1621)
        self.assertEqual(info.numCalpha(modelIdx=1), 218)

        return

    def testGetInfo2(self):
        """"""

        pdbfile = os.path.join(self.testfiles_dir, "2UUI.pdb")

        info = pdb_edit.get_info(pdbfile)

        self.assertEqual(len(info.models), 1)

        m1 = info.models[0]
        self.assertEqual(m1.chains[0], 'A')
        self.assertEqual(m1.resSeqs[0], [i for i in range(-5, 150)])
        self.assertEqual(m1.sequences[0],
                         'MHHHHHHKDEVALLAAVTLLGVLLQAYFSLQVISARRAFRVSPPLTTGPPEFERVYRAQVNCSEYFPLFLATLWVAGIFFHEGAAALCGLVYL'
                         'FARLRYFQGYARSAQLRLAPLYASARALWLLVALAALGLLAHFLPAALRAALLGRLRTLLPW')
        self.assertEqual(m1.caMask[0], [False] * 154 + [True])
        self.assertEqual(m1.bbMask[0],
                         [False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, True, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, True, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, True, False,
                          False, False, False, False, True, False, False, False, False, False, True, False, False,
                          False, False, False, False, False, False, False, False, False, False, True, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, False, False, False, False, False, False, False, False, False, False,
                          True, False, False, False, False, False, False, False, False, False, False, False, False,
                          False, False, False, True, False, False, False, False, True, True, True, True])

        self.assertEqual(info.numAtoms(modelIdx=0), 1263)

        return

    def testCheckPdbs(self):
        logging.basicConfig()
        logging.getLogger().setLevel(logging.DEBUG)

        pdbs = glob.glob(os.path.join(self.testfiles_dir, "models", "*.pdb"))
        self.assertTrue(pdb_edit.check_pdbs(pdbs))

        self.assertFalse(pdb_edit.check_pdbs(pdbs, single=True, sequence="AABBCC"))

        pdbs += [os.path.join(self.testfiles_dir, "1GU8.pdb")]
        self.assertFalse(pdb_edit.check_pdbs(pdbs, single=True, sequence="AABBCC"))

        return

    def testSelectResidues(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "testSelectResidues1.pdb"
        to_delete = [5, 10, 15, 20]

        b4 = set(pdb_edit.resseq(pdbin)['A'])

        pdb_edit.select_residues(pdbin=pdbin, pdbout=pdbout, delete=to_delete)

        after = set(pdb_edit.resseq(pdbout)['A'])
        self.assertEqual(after, b4.difference(set(to_delete)))

        os.unlink(pdbout)

        return

    def testSequence1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        ref = {'A': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
               'B': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
               'C': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY'}
        s = pdb_edit.sequence(pdbin)
        self.assertEqual(ref, s, "Bad _sequecne: {0}".format(s))
        return

    def XtestSplit(self):
        pdbin = os.path.join(self.testfiles_dir, "1GU8.pdb")
        pdb_edit.Xsplit(pdbin)
        # os.unlink(pdbout)
        return

    def testStripHetatm(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_het.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        pdb_edit._strip(hierachy, hetatm=True, hydrogen=False)
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('HETATM')])
        self.assertFalse(got, "Found HETATMS")
        os.unlink(pdbout)
        return

    def testStripHydrogen(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_H.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        pdb_edit._strip(hierachy, hetatm=False, hydrogen=True)
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('ATOM') and l[13] == 'H'])
        self.assertFalse(got, "Found Hydrogens")
        os.unlink(pdbout)
        return

    def testStripAtomTypes(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_types.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        pdb_edit._strip(hierachy, hetatm=False, hydrogen=False, atom_types=['CB'])
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('ATOM') and l[12:15].strip() == 'CB'])
        self.assertFalse(got, "Found Atom Types")
        os.unlink(pdbout)
        return

    def testReliableSidechains(self):
        pdbin = os.path.join(self.testfiles_dir, "1GU8.pdb")
        pdbout = "std.pdb"

        pdb_edit.reliable_sidechains(pdbin, pdbout)

        # Check it's valid
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)

        # Get list of all the residue names in chain 1
        resnames = [g.unique_resnames()[0] for g in pdb_obj.hierarchy.models()[0].chains()[0].residue_groups()]
        ref = ['VAL', 'GLY', 'LEU', 'THR', 'THR', 'LEU', 'PHE', 'TRP', 'LEU', 'GLY', 'ALA', 'ILE', 'GLY', 'MET',
               'LEU', 'VAL', 'GLY', 'THR', 'LEU', 'ALA', 'PHE', 'ALA', 'TRP', 'ALA', 'GLY', 'ARG', 'ASP', 'ALA',
               'GLY', 'SER', 'GLY', 'GLU', 'ARG', 'ARG', 'TYR', 'TYR', 'VAL', 'THR', 'LEU', 'VAL', 'GLY', 'ILE',
               'SER', 'GLY', 'ILE', 'ALA', 'ALA', 'VAL', 'ALA', 'TYR', 'VAL', 'VAL', 'MET', 'ALA', 'LEU', 'GLY',
               'VAL', 'GLY', 'TRP', 'VAL', 'PRO', 'VAL', 'ALA', 'GLU', 'ARG', 'THR', 'VAL', 'PHE', 'ALA', 'PRO',
               'ARG', 'TYR', 'ILE', 'ASP', 'TRP', 'ILE', 'LEU', 'THR', 'THR', 'PRO', 'LEU', 'ILE', 'VAL', 'TYR',
               'PHE', 'LEU', 'GLY', 'LEU', 'LEU', 'ALA', 'GLY', 'LEU', 'ASP', 'SER', 'ARG', 'GLU', 'PHE', 'GLY',
               'ILE', 'VAL', 'ILE', 'THR', 'LEU', 'ASN', 'THR', 'VAL', 'VAL', 'MET', 'LEU', 'ALA', 'GLY', 'PHE',
               'ALA', 'GLY', 'ALA', 'MET', 'VAL', 'PRO', 'GLY', 'ILE', 'GLU', 'ARG', 'TYR', 'ALA', 'LEU', 'PHE',
               'GLY', 'MET', 'GLY', 'ALA', 'VAL', 'ALA', 'PHE', 'LEU', 'GLY', 'LEU', 'VAL', 'TYR', 'TYR', 'LEU',
               'VAL', 'GLY', 'PRO', 'MET', 'THR', 'GLU', 'SER', 'ALA', 'SER', 'GLN', 'ARG', 'SER', 'SER', 'GLY',
               'ILE', 'LYS', 'SER', 'LEU', 'TYR', 'VAL', 'ARG', 'LEU', 'ARG', 'ASN', 'LEU', 'THR', 'VAL', 'ILE',
               'LEU', 'TRP', 'ALA', 'ILE', 'TYR', 'PRO', 'PHE', 'ILE', 'TRP', 'LEU', 'LEU', 'GLY', 'PRO', 'PRO',
               'GLY', 'VAL', 'ALA', 'LEU', 'LEU', 'THR', 'PRO', 'THR', 'VAL', 'ASP', 'VAL', 'ALA', 'LEU', 'ILE',
               'VAL', 'TYR', 'LEU', 'ASP', 'LEU', 'VAL', 'THR', 'LYS', 'VAL', 'GLY', 'PHE', 'GLY', 'PHE', 'ILE',
               'ALA', 'LEU', 'ASP', 'ALA', 'ALA', 'ALA', 'THR', 'LEU']

        self.assertEqual(resnames, ref)

        pdb_edit.reliable_sidechains_cctbx(pdbin, pdbout)
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)
        self.assertEqual(resnames, ref)
        os.unlink(pdbout)

        return

    def test_translate(self):
        # Test a translation function on a pdb file

        #######################################################################
        # Test case 1
        #######################################################################
        ftranslate = [1, 2, -1]
        pdbin = os.path.join(self.testfiles_dir, "2UUI.pdb")
        pdbout = "tmp_2UUI.pdb"

        reference_data = ["ATOM      1  N   MET A  -5     212.925 283.416-201.620  1.00 60.47           N",
                          "ATOM      2  CA  MET A  -5     214.345 283.308-201.183  1.00 60.38           C",
                          "ATOM      3  C   MET A  -5     215.036 284.648-201.449  1.00 59.70           C",
                          "ATOM      4  O   MET A  -5     216.159 284.672-201.948  1.00 60.50           O",
                          "ATOM      5  CB  MET A  -5     215.067 282.140-201.936  1.00 59.91           C",
                          "ATOM      6  N   HIS A  -4     214.350 285.747-201.100  1.00 58.69           N",
                          "ATOM      7  CA  HIS A  -4     214.689 287.114-201.530  1.00 56.01           C",
                          "ATOM      8  C   HIS A  -4     214.726 287.089-203.061  1.00 54.60           C",
                          "ATOM      9  O   HIS A  -4     215.327 286.191-203.618  1.00 54.54           O",
                          "ATOM     10  CB  HIS A  -4     216.035 287.562-200.905  1.00 56.11           C"
                          ]

        pdb_edit.translate(pdbin, pdbout, ftranslate)

        count = 0
        data = []
        with open(pdbout) as f:
            for line in f:
                if count < 10:
                    # Make sure the water and anisou atoms are deleted
                    if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
                        data.append(line[0:78])
                        count += 1
                    else:
                        continue
        os.unlink(pdbout)
        self.assertEqual(data, reference_data)
        return

    def testXyzCoordinates(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        test_hierarchy = iotbx.pdb.pdb_input(file_name=pdbin).construct_hierarchy()
        xyz_lst = pdb_edit._xyz_coordinates(test_hierarchy)

        ref_data_start = [(0, [(25.199, 11.913, -9.25),
                               (25.201, 10.666, -9.372),
                               (26.454, 12.702, -9.001)]),
                          (1, [(24.076, 12.643, -9.179),
                               (22.806, 12.124, -9.698),
                               (22.170, 11.067, -8.799),
                               (22.404, 11.024, -7.580)]),
                          (2, [(21.377, 10.190, -9.397),
                               (20.675, 9.156, -8.637),
                               (21.614, 8.106, -7.996),
                               (21.337, 7.619, -6.898),
                               (19.625, 8.485, -9.531),
                               (18.637, 7.595, -8.790),
                               (17.652, 8.361, -7.951),
                               (17.724, 9.603, -7.887),
                               (16.786, 7.706, -7.365)])]

        for idx in xrange(len(ref_data_start)):  # Stuff that needs to be true
            self.assertEqual(ref_data_start[idx][0], xyz_lst[idx][0])
            self.assertSequenceEqual(ref_data_start[idx][1], xyz_lst[idx][1])
        nr_atoms = sum(len(i[1]) for i in xyz_lst)
        self.assertEqual(252, nr_atoms)
        self.assertEqual(35, len(xyz_lst))

    def testXyzCbCoordinates(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        test_hierarchy = iotbx.pdb.pdb_input(file_name=pdbin).construct_hierarchy()
        xyz_cb_lst = pdb_edit._xyz_cb_coordinates(test_hierarchy)

        ref_data_start = [(0, (float('inf'), float('inf'), float('inf'))),
                          (1, (22.806, 12.124, -9.698)),
                          (2, (19.625, 8.485, -9.531)),
                          (3, (24.783, 6.398, -9.051)),
                          (4, (25.599, 10.846, -6.036)),
                          (5, (20.430, 10.143, -4.644))]

        self.assertSequenceEqual(ref_data_start[1], xyz_cb_lst[1][:6])
        self.assertEqual(35, len(xyz_cb_lst))


if __name__ == "__main__":
    unittest.main(verbosity=2)
