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
from rdkit import DataStructs
import argparse, textwrap

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

def calculateSimilarity(reference, ligand, method):
	if   method == 1:
		return DataStructs.TanimotoSimilarity(reference,ligand)
	elif method == 2:
		return DataStructs.DiceSimilarity(reference,ligand)
	elif method == 3:
		return DataStructs.CosineSimilarity(reference,ligand)
	elif method == 4:
		return DataStructs.SokalSimilarity(reference,ligand)
	elif method == 5:
		return DataStructs.RusselSimilarity(reference, ligand)
	elif method == 6:
		return DataStructs.RogotGoldbergSimilarity(reference, ligand)
	elif method == 7:
		return DataStructs.AllBitSimilarity(reference, ligand)
	elif method == 8:
		return DataStructs.KulczynskiSimilarity(reference, ligand)
	elif method == 9:
		return DataStructs.McConnaugheySimilarity(reference, ligand)
	elif method == 10:
		return DataStructs.AsymmetricSimilarity(reference, ligand)
	elif method == 11:
		return DataStructs.BraunBlanquetSimilarity(reference, ligand)
	elif method == 12:
		return DataStructs.TverskySimilarity(reference,ligand, args.alpha, args.beta)

if __name__ == '__main__':

	description = 'ProLIF: Protein Ligand Interaction Fingerprints\nGenerates Interaction FingerPrints (IFP) and a similarity score for protein-ligand interactions'
	epilog = 'Mandatory arguments: --reference --ligand --protein --residues\n1 single molecule per file.'
	# Argparse
	parser = argparse.ArgumentParser(description=description, epilog=epilog, 
		formatter_class=argparse.RawTextHelpFormatter)
	group_input = parser.add_argument_group('INPUT arguments')
	group_input.add_argument("-r", "--reference", metavar='fileName', type=str, required=True, 
		help="Path to your reference ligand. Type: SDF, MOL2, PDB")
	group_input.add_argument("-l", "--ligand", metavar='fileName', type=str, nargs='+', required=True, 
		help="Path to your ligand(s). Type: SDF, MOL2, PDB")
	group_input.add_argument("-p", "--protein", metavar='fileName', type=str, required=True, 
		help="Path to your protein or bindingSite. Type: PDB")
	group_input.add_argument("--residues", metavar='int', type=int, nargs='+', required=True, 
		help="Residues chosen for the interactions. Example: 11 32 34 171 174")
	group_input.add_argument("--chain", metavar="letter", type=str, default=None,
		help="Restricts the IFP to the residues present in specified chain.")
	group_output = parser.add_argument_group('OUTPUT arguments')
	group_output.add_argument("-o", "--output", metavar='filename', required=False, type=str, 
		help="Path to the output CSV file")
	group_output.add_argument("-v", "--verbose", action="store_true", help="Increase terminal output verbosity")
	group_args = parser.add_argument_group('Other arguments')
	group_args.add_argument("--similarity", choices=range(1,13), default=1, type=int, 
		help=textwrap.dedent("""Similarity score between molecule A and B :
Let 'a' and 'b' be the number of bits activated in molecules A and B,
'c' the number of activated bits in common, and 'd' the number of inactivated bits in common.
    1) Tanimoto      : c/(a+b-c). Used by default
    2) Dice          : 2c/(a+b)
    3) Cosine        : c/sqrt(a*b)
    4) Sokal         : c/(2a+2b-3c)
    5) Russel        : c/(a+b+d-c)
    6) RogotGoldberg : (c/(a+b))+(d/(a+b+2d-2c))
    7) AllBit        : (a - (a XOR b))/a
    8) Kulczynski    : c*(a+b)/(2*a*b)
    9) McConnaughey  : (c*(a+b)-(a*b))/(a*b)
   10) Asymmetric    : c/min(a,b)
   11) BraunBlanquet : c/max(a,b)
   12) Tversky       : c/(alpha*(a-c)+beta*(b-c)+c)
"""))
	group_args.add_argument("--alpha", metavar="int", type=float, choices=[Range(0.0, 1.0)], 
		help="Alpha parameter for Tversky")
	group_args.add_argument("--beta", metavar="int", type=float, choices=[Range(0.0, 1.0)], 
		help="Beta parameter for Tversky")

	# Parse arguments from command line
	args = parser.parse_args()
	# Read files
	reference = Ligand(args.reference)
	protein = Protein(args.protein, chain=args.chain, residueList=args.residues)
	# Print residues on terminal:
	for residue in protein.residues:
		print('{: >8s}'.format(protein.residues[residue].resName), end='')
	print()
	# Generate the IFP between the reference ligand and the protein
	reference.generateIFP(protein)
	print(reference.IFP)
	# Loop over ligands:
	ligandList = []
	for lig in args.ligand:
		ligand = Ligand(lig)
		# Generate the IFP between a ligand and the protein
		ligand.generateIFP(protein)
		# Calculate similarity
		S = calculateSimilarity(reference.IFPvector, ligand.IFPvector, args.similarity)
		ligand.setSimilarity(S)
		print(ligand.IFP, '{:.4f}'.format(ligand.S))
		ligandList.append(ligand)
	# Output
	if args.output:
		with open(args.output, 'w') as file:
			file.write('File,SimilarityScore')
			for residue in protein.residues:
				file.write(',{}'.format(protein.residues[residue].resName))
			file.write('\n')
			CSIFP = ','.join(reference.IFP[i:i+8] for i in range(0, len(reference.IFP), 8))
			file.write('{},1.0000,{}\n'.format(args.reference, CSIFP))
			for ligand in ligandList:
				CSIFP = ','.join(ligand.IFP[i:i+8] for i in range(0, len(ligand.IFP), 8))
				file.write('{},{:.4f},{}\n'.format(ligand.file, ligand.S, CSIFP))