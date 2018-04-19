import unittest
from os import path
from rdkit import DataStructs
from prolif.ligand import *
from prolif.fingerprint import Fingerprint

class TestLigand(unittest.TestCase):
    """Test the ligand.py module"""

    def setUp(self):
        self.here = path.abspath(path.dirname(__file__))
        self.lig  = Ligand(path.join(self.here,'files','dummy.mol2'))
        self.fp   = Fingerprint()

    def test_init_notmol2(self):
        with self.assertRaises(ValueError):
            lig = Ligand('ligand_test.py')

    def test_init_centroid(self):
        self.assertListEqual(self.lig.centroid, [-8.611,15.060,27.954])

    def test_repr(self):
        self.assertEqual(str(self.lig), path.join(self.here,'files','dummy.mol2'))

    def test_setIFP(self):
        self.lig.setIFP('101', 'vector')
        self.assertEqual(self.lig.IFP, '101')
        self.assertEqual(self.lig.IFPvector, 'vector')

    def test_setSimilarity(self):
        self.lig.setSimilarity(1.0)
        self.assertEqual(self.lig.score, 1.0)

    def test_getSimilarity(self):
        IFPvector = DataStructs.ExplicitBitVect(len(self.fp.interactions)*1)
        self.lig.setIFP('0', IFPvector)
        self.assertEqual(self.lig.getSimilarity(self.lig),1.0)

if __name__ == '__main__':
    unittest.main()
