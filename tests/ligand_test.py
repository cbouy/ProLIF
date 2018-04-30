import unittest
from os import path
from rdkit import DataStructs
from prolif.ligand import Ligand
from prolif.fingerprint import Fingerprint

class TestLigand(unittest.TestCase):
    """Test the ligand.py module"""

    def setUp(self):
        self.examples = path.abspath(path.join(path.dirname(__file__), '..', 'examples'))
        self.lig  = Ligand(path.join(self.examples,'ligand.mol2'))
        self.fp   = Fingerprint()


    def test_init_notmol2(self):
        with self.assertRaises(ValueError):
            lig = Ligand('ligand_test.py')

    def test_init_centroid(self):
        self.assertListEqual(self.lig.centroid, [-0.6397851063829787, 26.40136382978724, 10.785044680851064])


    def test_repr(self):
        self.assertEqual(str(self.lig), path.join(self.examples,'ligand.mol2'))


    def test_setIFP(self):
        self.lig.setIFP('101', 'vector')
        self.assertEqual(self.lig.IFP, '101')
        self.assertEqual(self.lig.IFPvector, 'vector')


    def test_setSimilarity(self):
        self.lig.setSimilarity(1.0)
        self.assertEqual(self.lig.score, 1.0)


    def test_getSimilarity_tanimoto(self):
        IFPvector = DataStructs.ExplicitBitVect(len(self.fp.interactions)*1)
        IFPvector.SetBit(0)
        self.lig.setIFP('1', IFPvector)
        self.assertEqual(self.lig.getSimilarity(self.lig, method='tanimoto'),1.0)

    def test_getSimilarity_dice(self):
        IFPvector = DataStructs.ExplicitBitVect(len(self.fp.interactions)*1)
        IFPvector.SetBit(0)
        self.lig.setIFP('1', IFPvector)
        self.assertEqual(self.lig.getSimilarity(self.lig, method='dice'),1.0)

    def test_getSimilarity_tversky(self):
        IFPvector = DataStructs.ExplicitBitVect(len(self.fp.interactions)*1)
        IFPvector.SetBit(0)
        self.lig.setIFP('1', IFPvector)
        self.assertEqual(self.lig.getSimilarity(self.lig, method='tversky', alpha=1, beta=0),1.0)

if __name__ == '__main__':
    unittest.main()
