#!/usr/bin/python3

"""
   Copyright 2017 Cédric BOUYSSET

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from src.ligand import *
from src.protein import *
from src.fingerprint import *
from src.utils import Range
import argparse, textwrap, os, sys

if __name__ == '__main__':

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
    group_input.add_argument("--json", metavar='fileName', type=str, default=os.path.join(os.path.dirname(sys.argv[0]),'prolif.json'),
        help="Path to a custom parameters file. Default: prolif.json")

    group_output = parser.add_argument_group('OUTPUT arguments')
    group_output.add_argument("-o", "--output", metavar='filename', type=str,
        help="Path to the output CSV file")
    group_output.add_argument("-v", "--verbose", action="store_true", help="Increase terminal output verbosity")

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
    group_args.add_argument("--similarity", choices=['tanimoto', 'dice', 'tversky'], default='tanimoto',
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
    # Read files
    fingerprint = Fingerprint(args.json, args.interactions)
    reference   = Ligand(args.reference)
    protein     = Protein(args.protein, reference, cutoff=args.cutoff, residueList=args.residues)
    length      = len(args.interactions)
    # Print residues on terminal:
    for residue in protein.residues:
        print('{resid: >{length}s}'.format(resid=protein.residues[residue].resid, length=length), end='')
    print()
    # Generate the IFP between the reference ligand and the protein
    fingerprint.generateIFP(reference, protein)
    print(reference.IFP)
    # Loop over ligands:
    ligandList = []
    for lig in args.ligand:
        ligand = Ligand(lig)
        # Generate the IFP between a ligand and the protein
        fingerprint.generateIFP(ligand, protein)
        # Calculate similarity
        S = ligand.getSimilarity(reference, args.similarity, args.alpha, args.beta)
        ligand.setSimilarity(S)
        print(ligand.IFP, '{:.4f}'.format(ligand.S))
        ligandList.append(ligand)
    # Output
    if args.output:
        with open(args.output, 'w') as file:
            file.write('File,SimilarityScore')
            for residue in protein.residues:
                file.write(',{}'.format(protein.residues[residue].resid))
            file.write('\n')
            CSIFP = ','.join(reference.IFP[i:i+length] for i in range(0, len(reference.IFP), length))
            file.write('{},1.0000,{}\n'.format(args.reference, CSIFP))
            for ligand in ligandList:
                CSIFP = ','.join(ligand.IFP[i:i+length] for i in range(0, len(ligand.IFP), length))
                file.write('{},{:.4f},{}\n'.format(ligand.file, ligand.S, CSIFP))
