import unittest
from os import path
import subprocess
from prolif.version import __version__

class TestCLI(unittest.TestCase):
    """Test the command_line.py module. Needs prolif to be installed"""

    def test_version(self):
        process = subprocess.run(['prolif', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(process.stdout.decode('utf-8'), 'ProLIF {}\n'.format(__version__))

if __name__ == '__main__':
    unittest.main()
