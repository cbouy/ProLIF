import unittest
from os import path
from prolif.protein import Protein
from prolif.ligand import Ligand

class TestProtein(unittest.TestCase):
    """Test the protein.py module"""

    def setUp(self):
        self.examples = path.abspath(path.join(path.dirname(__file__), '..', 'examples'))
        self.reference = Ligand(path.join(self.examples,'ligand.mol2'))
        self.protein = Protein(
            inputFile=path.join(self.examples,'protein.mol2'),
            reference=self.reference,
            residueList=['LEU83'])


    def test_repr(self):
        self.assertEqual(str(self.protein), path.join(self.examples,'protein.mol2'))

    def test_init_notmol2(self):
        with self.assertRaises(ValueError):
            Protein('not a mol2 file')

    def test_init_detectCloseResidues(self):
        residues = self.protein.detectCloseResidues(self.reference)
        self.assertListEqual(residues, ['LEU83'])

    def test_init_detectCloseResidues_shortcutoff(self):
        residues = self.protein.detectCloseResidues(self.reference, cutoff=1.0)
        self.assertListEqual(residues, [])

    def test_cleanResidues(self):
        self.protein.residues = {'A':[],'B':[],'C':[]}
        self.protein.residueList = ['B']
        self.protein.cleanResidues()
        self.assertDictEqual(self.protein.residues, {'B':[]})


if __name__ == '__main__':
    unittest.main()
