import unittest
from prolif.residue import *

class TestResidue(unittest.TestCase):
    """Test the residue.py module"""

    def setUp(self):
        self.atoms = [
            {
                'resname':'ALA',
                'resid':'ALA1',
                'coordinates': [0,1,2]
            },
            {
                'resname':'ALA',
                'resid':'ALA1',
                'coordinates': [2,3,4]
            },
        ]
        self.residue  = Residue(self.atoms)


    def test_init_atoms(self):
        self.assertListEqual(self.atoms, self.residue.atoms)

    def test_init_resname(self):
        self.assertEqual(self.residue.resname, 'ALA')

    def test_init_resid(self):
        self.assertEqual(self.residue.resid, 'ALA1')

    def test_init_centroid(self):
        self.assertListEqual(self.residue.centroid, [1,2,3])


    def test_repr(self):
        self.assertEqual(str(self.residue), 'ALA1')

if __name__ == '__main__':
    unittest.main()
