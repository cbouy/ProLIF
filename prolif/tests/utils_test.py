import unittest
from os import path
from prolif.utils import *

class TestUtils(unittest.TestCase):
    """Test the utils.py module"""

    def test_euclidianDistance_0_0_0__1_0_0(self):
        self.assertEqual(euclidianDistance([0,0,0],[1,0,0]), 1)

    def test_euclidianDistance_0__0(self):
        self.assertNotEqual(euclidianDistance([0],[0]), 1)

    def test_euclidianDistance_a__0(self):
        with self.assertRaises(TypeError):
            euclidianDistance(['a'], [0])

    def test_getCentroid_empty_atoms(self):
        atoms = [
            {'coordinates':[]},
        ]
        with self.assertRaises(IndexError):
            getCentroid(atoms)

    def test_getCentroid_not3d(self):
        atoms = [
            {'coordinates':[0,0]},
            {'coordinates':[1,1]},
        ]
        with self.assertRaises(IndexError):
            getCentroid(atoms)

    def test_getCentroid_no_atoms(self):
        atoms = []
        with self.assertRaises(ZeroDivisionError):
            getCentroid(atoms)

    def test_getCentroid_0_0_0(self):
        atoms = [
            {'coordinates':[0,0,0]},
        ]
        self.assertEqual(getCentroid(atoms), [0,0,0])

    def test_getCentroid_0_0_0__2_2_2(self):
        atoms = [
            {'coordinates':[0,0,0]},
            {'coordinates':[2,2,2]},
        ]
        self.assertEqual(getCentroid(atoms), [1,1,1])

    def test_get_mol_from_mol2_nolines(self):
        self.assertEqual(get_mol_from_mol2(0,0,['0 1 2 3 4 5 6 7'],False), [])

    def test_get_mol_from_mol2_nb_atoms_not_int(self):
        with self.assertRaises(ValueError):
            get_mol_from_mol2(0,0,['foo 1 2 3 4 5 6 7'],False)

    def test_get_mol_from_mol2_hydrogen_ignored(self):
        self.assertEqual(get_mol_from_mol2(0,0,['1 1 2 3 4 H 6 7'],True), [])

    def test_get_mol_from_mol2_hydrogen_notignored(self):
        self.assertNotEqual(get_mol_from_mol2(0,0,['1 1 2 3 4 H 6 7'],False), [])

    def test_get_mol_from_mol2_not_hydrogen_notignored(self):
        self.assertNotEqual(get_mol_from_mol2(0,0,['1 1 2 3 4 C 6 7'],True), [])

    def test_get_mol_from_mol2_not_enough_fields(self):
        with self.assertRaises(IndexError):
            get_mol_from_mol2(0,0,['1 1 2 3 4'],False)

    def test_mol2_reader_nofile(self):
        with self.assertRaises(FileNotFoundError):
            mol2_reader('not a mol2 file', False)

    def test_mol2_reader_dummyfile(self):
        here = path.abspath(path.dirname(__file__))
        self.assertEqual(mol2_reader(path.join(here,'files','dummy.mol2'),False), [{'chain':'1','residue':'XK2263'}])

if __name__ == '__main__':
    unittest.main()