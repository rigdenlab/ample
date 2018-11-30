
import pickle
import os
import unittest
from ample.util import benchmark_util

@unittest.skip("not a regular test")
class Test(unittest.TestCase):
    def test_benchmark(self):
        pklfile="/home/jmht/ample-dev1/examples/toxd-example/ROSETTA_MR_0/resultsd.pkl"
        with open(pklfile) as f: d=pickle.load(f)
        bd="/home/jmht/ample-dev1/python/foo"
        if not os.path.isdir(bd): os.mkdir(bd)
        d['benchmark_dir']=bd
        benchmark_util.analyse(d)
