"""Test functions for util.ample_util"""
import unittest
from ample.util import reference_util


class Test(unittest.TestCase):

    def test_construct_references(self):
        #import argparse
        from ample.util import config_util, argparse_util
        options = config_util.AMPLEConfigOptions()
        argso = argparse_util.process_command_line(args=['-mtz', 'foo',
                                                         '-fasta', 'bar'])
        options.populate(argso)
        references = reference_util.construct_references(options.d, write_file=False)
        ref_references = '* Zhang et al. (2004). SPICKER: A clustering approach to identify near-native protein folds. Journal of Computational Chemistry 25(6), 865-871. [doi:10.1002/jcc.20011]'
        self.assertEqual(references, ref_references)
        
        options.d['nmr_model_in'] = 'foo'
        options.d['transmembrane'] = True
        options.d['use_scwrl'] = True
        options.d['do_mr'] = True
        options.d['refine_rebuild_arpwarp']  = True
        options.d['shelxe_rebuild_buccaneer'] = True
        options.d['use_shelxe'] = True
        options.d['mrbump_programs'] = ['molrep', 'phaser']
        ref_references = '* Bibby et al. (2012). AMPLE: A cluster-and-truncate approach to solve the crystal structures of small proteins using rapidly computed ab initio models. Acta Crystallogr. Sect. D Biol. Crystallogr. 68(12), 1622-1631. [doi:10.1107/S0907444912039194]\n\n* Bibby et al. (2013). Application of the AMPLE cluster-and-truncate approach to NMR structures for molecular replacement. Acta Crystallogr. Sect. D Biol. Crystallogr. 69(11), 2194-2201. [doi:10.1107/S0907444913018453]\n\n* Keegan et al. (2015). Exploring the speed and performance of molecular replacement with AMPLE using QUARK ab initio protein models. Acta Crystallogr. Sect. D Biol. Crystallogr. 71(2), 338-343. [doi:10.1107/S1399004714025784]\n\n* Thomas et al. (2015). Routine phasing of coiled-coil protein crystal structures with AMPLE. IUCrJ 2(2), 198-206. [doi:10.1107/S2052252515002080]\n\n* Simkovic et al. (2016). Residue contacts predicted by evolutionary covariance extend the application of ab initio molecular replacement to larger and more challenging protein folds. IUCrJ 3(4), 259-270. [doi:10.1107/S2052252516008113]\n\n* Thomas et al. (2017). Approaches to ab initio molecular replacement of alpha-helical transmembrane proteins. Acta Crystallographica Section D 73(12), 985-996. [doi:10.1107/S2059798317016436]\n\n* Bradley et al. (2005). Toward High-Resolution de Novo Structure Prediction for Small Proteins. Science 309(5742), 1868-1871. [doi:10.1126/science.1113801]\n\n* Theobald et al. (2006). THESEUS: maximum likelihood superpositioning and analysis of macromolecular structures. Bioinformatics 22(17), 2171-2172. [doi:10.1093/bioinformatics/btl332]\n\n* Krivov et al. (2009). Improved prediction of protein side-chain conformations with SCWRL4. Proteins: Struct., Funct., Bioinf. 77(4), 778-795. [doi:10.1002/prot.22488]\n\n* Winn et al. (2011). Overview of the CCP4 suite and current developments. Acta Crystallographica Section D 67(4), 235-242. [doi:10.1107/S0907444910045749]\n\n* Krissinel et al. (2012). Enhanced fold recognition using efficient short fragment clustering. Journal of molecular biochemistry 1(2), 76\xe2\x80\x9485. [doi:]'
        references = reference_util.construct_references(options.d, write_file=False)
        self.assertEqual(references, ref_references)

if __name__ == "__main__":
    unittest.main()
