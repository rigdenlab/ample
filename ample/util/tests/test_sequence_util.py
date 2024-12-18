"""Test functions for util.sequence_util"""

import os
import unittest
from ample import constants
from ample.util import sequence_util


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')

    def testSequence1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        ref = {
            'A': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
            'B': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
            'C': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
        }
        s = sequence_util.sequence(pdbin)
        self.assertEqual(ref, s, "Bad sequence: {0}".format(s))
        return

    def test_add(self):
        s1 = sequence_util.Sequence(pdb=os.path.join(self.testfiles_dir, '1GU8.pdb'))
        s2 = sequence_util.Sequence(fasta=os.path.join(self.testfiles_dir, '2uui.fasta'))
        s1 += s2

        self.assertTrue(len(s1.sequences), 2)
        self.assertTrue(len(s1.resseqs), 2)
        self.assertTrue(len(s1.headers), 2)
        self.assertTrue(len(s1.pdbs), 2)
        self.assertTrue(len(s1.chains), 2)
        self.assertTrue(len(s1.fasta_files), 2)

    def test_addPdb_data(self):
        fasta1 = os.path.join(self.testfiles_dir, '1ujb_2a6pA_3c7tA.afasta')
        pdbin1 = os.path.join(self.ample_share, 'examples', 'homologs', 'input', '1ujbA.pdb')
        pdbin2 = os.path.join(self.ample_share, 'examples', 'homologs', 'input', '2a6pA.pdb')
        pdbin3 = os.path.join(self.ample_share, 'examples', 'homologs', 'input', '3c7tA.pdb')
        s1 = sequence_util.Sequence(fasta=fasta1)
        s1.add_pdb_data(pdbin1)
        s1.add_pdb_data(pdbin2)
        s1.add_pdb_data(pdbin3)

        self.assertEqual(s1.pdbs[0], os.path.basename(pdbin1))
        self.assertEqual(s1.chains[0], 'A')
        self.assertEqual(s1.pdbs[1], os.path.basename(pdbin2))
        self.assertEqual(s1.chains[1], 'A')
        self.assertEqual(s1.pdbs[2], os.path.basename(pdbin3))
        self.assertEqual(s1.chains[2], 'A')

        p1r = [
            None,
            None,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            None,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            55,
            56,
            57,
            58,
            59,
            60,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            69,
            None,
            None,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            83,
            None,
            None,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            None,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            105,
            106,
            107,
            108,
            109,
            110,
            111,
            112,
            113,
            114,
            115,
            116,
            117,
            118,
            119,
            120,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            121,
            122,
            123,
            None,
            None,
            None,
            124,
            125,
            None,
            None,
            126,
            None,
            None,
            127,
            128,
            129,
            130,
            131,
            132,
            133,
            134,
            135,
            136,
            137,
            138,
            None,
            None,
            139,
            140,
            141,
            142,
            143,
            144,
            145,
            146,
            147,
            148,
            149,
            None,
            None,
            150,
            151,
            None,
            152,
            153,
            154,
            155,
            156,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ]

        p2r = [
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            23,
            24,
            25,
            26,
            27,
            28,
            None,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            None,
            51,
            52,
            53,
            54,
            55,
            56,
            57,
            58,
            59,
            60,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            69,
            70,
            71,
            72,
            73,
            None,
            None,
            None,
            None,
            74,
            None,
            None,
            75,
            None,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            105,
            106,
            107,
            108,
            109,
            110,
            111,
            112,
            113,
            114,
            115,
            116,
            117,
            118,
            119,
            120,
            121,
            122,
            123,
            124,
            125,
            126,
            127,
            128,
            129,
            130,
            131,
            132,
            None,
            133,
            134,
            135,
            136,
            137,
            None,
            None,
            138,
            139,
            140,
            141,
            142,
            143,
            144,
            145,
            146,
            147,
            148,
            149,
            150,
            151,
            152,
            153,
            154,
            155,
            156,
            157,
            158,
            159,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            160,
            161,
            162,
            163,
            None,
            164,
            165,
            166,
            167,
            168,
            169,
            None,
            None,
            170,
            171,
            172,
            173,
            174,
            175,
            176,
            177,
            178,
            179,
            180,
            181,
            182,
            183,
            184,
            None,
            None,
            None,
            185,
            186,
            187,
            188,
            189,
            190,
            191,
            192,
            None,
            193,
            194,
            195,
            196,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ]

        p3r = [
            None,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            105,
            106,
            107,
            108,
            109,
            110,
            111,
            112,
            113,
            114,
            115,
            116,
            117,
            118,
            119,
            120,
            121,
            122,
            123,
            124,
            None,
            125,
            126,
            127,
            128,
            129,
            130,
            131,
            132,
            133,
            134,
            135,
            136,
            137,
            138,
            139,
            140,
            141,
            142,
            143,
            144,
            145,
            146,
            147,
            148,
            149,
            None,
            150,
            151,
            152,
            153,
            154,
            155,
            156,
            157,
            158,
            159,
            160,
            161,
            162,
            163,
            164,
            165,
            166,
            167,
            168,
            169,
            170,
            171,
            172,
            173,
            174,
            175,
            176,
            177,
            178,
            179,
            180,
            181,
            182,
            183,
            184,
            185,
            186,
            187,
            188,
            189,
            190,
            191,
            192,
            193,
            194,
            195,
            196,
            197,
            198,
            199,
            200,
            201,
            202,
            203,
            204,
            205,
            206,
            207,
            208,
            209,
            210,
            211,
            212,
            213,
            214,
            215,
            216,
            217,
            218,
            219,
            220,
            221,
            222,
            223,
            224,
            225,
            226,
            227,
            228,
            229,
            230,
            231,
            232,
            233,
            234,
            235,
            236,
            237,
            238,
            239,
            240,
            241,
            242,
            243,
            244,
            245,
            246,
            247,
            248,
            249,
            250,
            None,
            251,
            252,
            253,
            254,
            255,
            256,
            257,
            258,
            259,
            260,
            261,
            262,
            263,
            264,
            265,
            266,
            267,
            268,
            269,
            270,
            271,
            272,
            273,
            274,
            275,
            276,
            277,
            278,
            279,
            280,
            281,
            282,
            283,
            284,
            285,
            286,
            287,
            288,
            None,
            None,
            289,
            290,
            291,
            292,
            293,
            294,
            295,
            296,
            297,
            298,
            299,
            300,
            301,
            302,
            303,
            None,
            None,
            304,
            None,
            None,
            None,
            305,
            306,
            307,
            308,
            309,
            310,
            311,
            312,
            313,
            314,
            315,
            None,
            316,
            317,
            318,
            319,
            320,
            321,
            322,
            323,
            324,
            325,
            326,
            327,
            328,
            329,
            330,
        ]

        self.assertEqual(s1.resseqs[0], p1r)
        self.assertEqual(s1.resseqs[1], p2r)
        self.assertEqual(s1.resseqs[2], p3r)

    def test_canonicalise_1(self):
        fp = sequence_util.Sequence()
        fp.sequences = ["YFLVKGMGVSDPDAKKFYAITTLVYAIAFTMYLSMLLGYGLTMVP"]
        fp.canonicalise()
        self.assertTrue(True)

    def test_canonicalise_2(self):
        fp = sequence_util.Sequence()
        fp.sequences = ["YFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVP"]
        with self.assertRaises(RuntimeError):
            fp.canonicalise()

    def test_canonicalise_3(self):
        fp = sequence_util.Sequence()
        fp.sequences = ["YFLVKGMG VSDPDAKKFY AITTLVXA IAFTMYLS MLLGYGLTMVP"]
        with self.assertRaises(RuntimeError):
            fp.canonicalise()

    def test_canonicalise_4(self):
        fp = sequence_util.Sequence()
        fp.sequences = ["YFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVP*"]
        with self.assertRaises(RuntimeError):
            fp.canonicalise()

    def test_OK(self):
        infasta = ">3HAP:A|PDBID|CHAIN|SEQUENCE" + os.linesep
        infasta += "QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAI" + os.linesep
        infasta += "TTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLV" + os.linesep
        infasta += "DADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAA" + os.linesep
        infasta += "MLYILYVLFFGFTSKAESMRPEVASTFKVL" + os.linesep
        infasta += "RNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILL" + os.linesep
        infasta += "RSRAIFGEAEAPEPSAGDGAAATSD"

        fp = sequence_util.Sequence()
        fp._parse_fasta(infasta.split(os.linesep))

        outfasta = ">3HAP:A|PDBID|CHAIN|SEQUENCE" + os.linesep
        outfasta += "QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYW" + os.linesep
        outfasta += "ARYADWLFTTPLLLLDLALLVDADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKA" + os.linesep
        outfasta += "ESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA" + os.linesep
        outfasta += "GDGAAATSD" + os.linesep
        outfasta += os.linesep

        self.assertEqual(outfasta, "".join(fp.fasta_str()))
        self.assertEqual(fp.length(), 249)

    def test_from_pdb(self):
        s1 = sequence_util.Sequence(pdb=os.path.join(self.testfiles_dir, '4DZN.pdb'))
        self.assertEqual(s1.name, '4DZN')
        self.assertEqual(s1.pdbs, ['4DZN.pdb', '4DZN.pdb', '4DZN.pdb'])
        self.assertEqual(s1.chains, ['A', 'B', 'C'])

        outfasta = ">From pdb: 4DZN.pdb chain=A length=31" + os.linesep
        outfasta += "GEIAALKQEIAALKKEIAALKEIAALKQGYY" + os.linesep
        outfasta += os.linesep
        outfasta += ">From pdb: 4DZN.pdb chain=B length=31" + os.linesep
        outfasta += "GEIAALKQEIAALKKEIAALKEIAALKQGYY" + os.linesep
        outfasta += os.linesep
        outfasta += ">From pdb: 4DZN.pdb chain=C length=31" + os.linesep
        outfasta += "GEIAALKQEIAALKKEIAALKEIAALKQGYY" + os.linesep
        outfasta += os.linesep

        self.assertEqual(outfasta, "".join(s1.fasta_str()))

    def test_resseq(self):
        pdbin = os.path.join(self.testfiles_dir, '1D7M.pdb')
        s1 = sequence_util.Sequence(pdb=pdbin)
        self.assertTrue(len(s1.sequences), 2)
        self.assertTrue(len(s1.headers), 2)
        self.assertTrue(len(s1.pdbs), 2)
        self.assertEqual(s1.pdbs[0], os.path.basename(pdbin))
        self.assertTrue(s1.resseqs[0][-1], 343)

    def test_align_file(self):
        pdbin1 = os.path.join(self.testfiles_dir, '1D7M.pdb')
        pdbin2 = os.path.join(self.testfiles_dir, '1GU8.pdb')
        pdbin3 = os.path.join(self.testfiles_dir, '2UUI.pdb')
        s1 = sequence_util.Sequence(pdb=pdbin1)
        s1 += sequence_util.Sequence(pdb=pdbin2)
        s1 += sequence_util.Sequence(pdb=pdbin3)

        ref = ">1D7M.pdb" + os.linesep
        ref += "EMANRLAGLENSLESEKVSREQLIKQKDQLNSLLASLESEGAEREKRLRELEAKLDETLKNLELEKLARMELEARLAKTE" + os.linesep
        ref += "KDRAILELKLAEAIDEKSKLE" + os.linesep
        ref += os.linesep
        ref += ">1D7M.pdb" + os.linesep
        ref += "EMANRLAGLENSLESEKVSREQLIKQKDQLNSLLASLESEGAEREKRLRELEAKLDETLKNLELEKLARMELEARLAKTE" + os.linesep
        ref += "KDRAILELKLAEAIDEKSKLE" + os.linesep
        ref += os.linesep
        ref += ">1GU8.pdb" + os.linesep
        ref += "VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTP" + os.linesep
        ref += "LIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYV" + os.linesep
        ref += "RLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATL" + os.linesep
        ref += os.linesep
        ref += ">2UUI.pdb" + os.linesep
        ref += "MHHHHHHKDEVALLAAVTLLGVLLQAYFSLQVISARRAFRVSPPLTTGPPEFERVYRAQVNCSEYFPLFLATLWVAGIFF" + os.linesep
        ref += "HEGAAALCGLVYLFARLRYFQGYARSAQLRLAPLYASARALWLLVALAALGLLAHFLPAALRAALLGRLRTLLPWA" + os.linesep
        ref += os.linesep

        self.assertEqual(s1.fasta_str(pdbname=True), ref)

    def test__parse_fasta_1(self):
        fasta = [">foo"]
        fasta += ["AAAAAAA"]
        s = sequence_util.Sequence()
        s._parse_fasta(fasta)
        self.assertListEqual(s.headers, [">foo"])
        self.assertListEqual(s.sequences, ["AAAAAAA"])

    def test__parse_fasta_2(self):
        fasta = [">foo"]
        fasta += ["AAAAA AA"]
        s = sequence_util.Sequence()
        s._parse_fasta(fasta)
        self.assertListEqual(s.headers, [">foo"])
        self.assertListEqual(s.sequences, ["AAAAAAA"])

    def test__parse_fasta_3(self):
        fasta = [">foo"]
        fasta += ["AAAAAAA*"]
        s = sequence_util.Sequence()
        s._parse_fasta(fasta)
        self.assertListEqual(s.headers, [">foo"])
        self.assertListEqual(s.sequences, ["AAAAAAA"])

    def test__parse_fasta_4(self):
        fasta = [">foo"]
        fasta += ["AAAA AAA"]
        fasta += [">bar"]
        fasta += ["CCCCCCC*"]
        s = sequence_util.Sequence()
        s._parse_fasta(fasta)
        self.assertListEqual(s.headers, [">foo", ">bar"])
        self.assertListEqual(s.sequences, ["AAAAAAA", "CCCCCCC"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
