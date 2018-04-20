import unittest
from os import path
from prolif.protein import *
from prolif.ligand import Ligand

class TestProtein(unittest.TestCase):
    """Test the protein.py module"""

    def setUp(self):
        self.here = path.abspath(path.dirname(__file__))
        self.reference = Ligand(path.join(self.here,'files','dummy.mol2'))
        self.protein = Protein(
            inputFile=path.join(self.here,'files','1hvr_rec.mol2'),
            reference=self.reference,
            residueList=['PRO1'])


    def test_repr(self):
        self.assertEqual(str(self.protein), path.join(self.here,'files','1hvr_rec.mol2'))

    def test_init_notmol2(self):
        with self.assertRaises(ValueError):
            Protein('not a mol2 file')

    def test_init_detectCloseResidues(self):
        residues = self.protein.detectCloseResidues(self.reference, cutoff=40.0)
        self.assertListEqual(residues, ['PRO1'])

    def test_init_detectCloseResidues_shortcutoff(self):
        residues = self.protein.detectCloseResidues(self.reference, cutoff=1.0)
        self.assertListEqual(residues, [])

    def test_cleanResidues(self):
        self.protein.residues = {'A':[],'B':[],'C':[]}
        self.protein.residueList = ['B']
        self.protein.cleanResidues()
        self.assertDictEqual(self.protein.residues, {'B':[]})

    def test_residuesFromMOL2File(self):
        self.protein.inputFile =  path.join(self.here,'files','dummy.mol2')
        self.protein.residuesFromMOL2File()
        self.assertEqual(self.protein.residues['XK2263'].resid, 'XK2263')

if __name__ == '__main__':
    unittest.main()
