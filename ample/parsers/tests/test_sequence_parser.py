
import unittest
from ample.parsers import sequence_parser

class Test(unittest.TestCase):
    
    def test_canonicalise(self):
        in_seq = "MAAAKDEVALLAAVTLLGVLL"
        seq = sequence_parser.canonicalise(in_seq)
        self.assertEqual(seq, "MAAAKDEVALLAAVTLLGVLL")
        
        in_seq = "MXXXKDEVALLAAVTLLGVLL"
        with self.assertRaises(RuntimeError):
            seq = sequence_parser.canonicalise(in_seq)
            
        in_seq = "MaaaKDEVALLAAVTLLGVLL"
        seq = sequence_parser.canonicalise(in_seq)
        self.assertEqual(seq, "MAAAKDEVALLAAVTLLGVLL")

if __name__ == "__main__":
    unittest.main()