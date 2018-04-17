import argparse, textwrap, os, sys
from . import prolif
from .version import __version__

class Range:
    """Class to raise an exception if the value is not between start and end.
    Used for alpha and beta in Tversky"""
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{} to {}'.format(self.start, self.end)

def cli():
    jsonpath = os.path.join(os.path.dirname(__file__),'prolif.json')
    description = 'ProLIF: Protein Ligand Interaction Fingerprints\nGenerates Interaction FingerPrints (IFP) and a similarity score for protein-ligand interactions'
    epilog = 'Mandatory arguments: --reference --ligand --protein\nMOL2 files only.'
    # Argparse
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter)

    group_input = parser.add_argument_group('INPUT arguments')
    group_input.add_argument("-r", "--reference", metavar='fileName', type=str, required=True,
        help="Path to your reference ligand.")
    group_input.add_argument("-l", "--ligand", metavar='fileName', type=str, nargs='+', required=True,
        help="Path to your ligand(s).")
    group_input.add_argument("-p", "--protein", metavar='fileName', type=str, required=True,
        help="Path to your protein.")
    group_input.add_argument("--residues", type=str, nargs='+', default=None,
        help="Residues chosen for the interactions. Default: automatically detect residues within --cutoff of the reference ligand")
    group_input.add_argument("--cutoff", metavar='float', type=float, required=False, default=5.0,
        help="Cutoff for automatic residue detection. Default: 5.0 angströms")
    group_input.add_argument("--json", metavar='fileName', type=str, default=jsonpath,
        help="Path to a custom parameters file. Default: prolif.json")

    group_output = parser.add_argument_group('OUTPUT arguments')
    group_output.add_argument("-o", "--output", metavar='filename', type=str,
        help="Path to the output CSV file")
    group_output.add_argument("-v", "--verbose", action="store_true", help="Increase terminal output verbosity")
    group_output.add_argument("--version", action="version",
        version='ProLIF {}'.format(__version__) ,help="Show version and exit")

    group_args = parser.add_argument_group('Other arguments')
    group_args.add_argument("--interactions", metavar="bit", nargs='+',
        choices=['HBdonor','HBacceptor','cation','anion','FaceToFace','FaceToEdge',
            'pi-cation','hydrophobic','XBdonor','metal'],
        default=['HBdonor','HBacceptor','cation','anion','FaceToFace','FaceToEdge','pi-cation','hydrophobic'],
        help=textwrap.dedent("""List of interactions used to build the fingerprint.
    -) hydrogen bond: HBdonor, HBacceptor
    -) halogen bond:  XBdonor
    -) ionic: cation, anion
    -) pi-stacking: FaceToFace, FaceToEdge
    -) hydrophobic
    -) pi-cation
    -) metal
    Default: HBdonor HBacceptor cation anion FaceToFace FaceToEdge pi-cation hydrophobic"""))
    group_args.add_argument("--score", choices=['tanimoto', 'dice', 'tversky'], default='tanimoto',
        help=textwrap.dedent("""Similarity score between molecule A and B :
    Let 'a' and 'b' be the number of bits activated in molecules A and B, and 'c' the number of activated bits in common.
    -) tanimoto : c/(a+b-c). Used by default
    -) dice     : 2c/(a+b)
    -) tversky  : c/(alpha*(a-c)+beta*(b-c)+c)
    """))
    group_args.add_argument("--alpha", metavar="int", type=float, choices=[Range(0.0, 1.0)], default=0.7,
        help="Alpha parameter for Tversky. Default: 0.7")
    group_args.add_argument("--beta", metavar="int", type=float, choices=[Range(0.0, 1.0)], default=0.3,
        help="Beta parameter for Tversky. Default: 0.3")

    # Parse arguments from command line
    args = parser.parse_args()

    prolif.main(args)
