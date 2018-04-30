import unittest
from os import path
from prolif.utils import *

class TestUtils(unittest.TestCase):
    """Test the utils.py module"""

    def test_euclidianDistance_0_0_0__1_0_0(self):
        self.assertEqual(euclidianDistance([0,0,0],[1,0,0]), 1)

    def test_euclidianDistance_0__0(self):
        self.assertNotEqual(euclidianDistance([0],[0]), 1)

    def test_euclidianDistance_0_0__0(self):
        with self.assertRaises(IndexError):
            euclidianDistance([0,0], [0])


    def test_getCentroid_empty_coords(self):
        with self.assertRaises(IndexError):
            getCentroid([[]])

    def test_getCentroid_not3d(self):
        with self.assertRaises(IndexError):
            getCentroid([[0,0],])

    def test_getCentroid_0_0_0(self):
        self.assertListEqual(getCentroid([[0,0,0],]), [0,0,0])

    def test_getCentroid_0_0_0__2_2_2(self):
        coords = [
            [0,0,0],
            [2,2,2],
        ]
        self.assertListEqual(getCentroid(coords), [1,1,1])


    def test_mol2_reader_nofile(self):
        with self.assertRaises(FileNotFoundError):
            mol2_reader('not a mol2 file')


if __name__ == '__main__':
    unittest.main()
