import unittest
from os import path
from rdkit import Chem
from prolif.residue import *

class TestResidue(unittest.TestCase):
    """Test the residue.py module"""

    def setUp(self):
        self.mol = Chem.MolFromMol2File(path.join(path.abspath(path.dirname(__file__)), 'files', 'dummy.mol2'))
        self.mol.SetProp('resname', 'ALA1')
        self.residue = Residue(self.mol)


    def test_init_mol(self):
        self.assertEqual(self.mol, self.residue.mol)

    def test_init_resname(self):
        self.assertEqual(self.residue.resname, 'ALA1')


    def test_repr(self):
        self.assertEqual(str(self.residue), 'ALA1')

if __name__ == '__main__':
    unittest.main()
