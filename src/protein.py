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

from rdkit import Chem
import os.path
from src.residue import *

class Protein:
	"""Class for a protein"""
	
	def __init__(self, file, residueList, chain):
		"""Initialization of the protein, defined by a list of residues"""
		self.residueList = residueList
		self.chain = chain
		self.file = file
		fileExtension = os.path.splitext(file)[1]
		if fileExtension == '.pdb':
			self.residuesFromPDBFile()
		else:
			raise ValueError('{} files are not supported for the protein.'.format(fileExtension[1:].upper()))

	def residuesFromPDBFile(self):
		"""Read a PDB file and assign each line to an object of class Atom"""
		# If a specific chain has been set:
		if self.chain:
			# Use RDKit to keep only the asked chain
			chains = Chem.SplitMolByPDBChainId(Chem.MolFromPDBFile(self.file))
			newFilename = '{}.chain{}'.format(self.file, self.chain)
			Chem.MolToPDBFile(chains[self.chain], newFilename)
			self.file = newFilename
		# Read the file
		with open(self.file, 'r') as file:
			# initialization
			residuesDict = {}
			for residue in self.residueList:
				residuesDict[residue] = []
			# read each line and assign atoms of residues to a dictionary
			for line in file.readlines():
				if line[0:4] == 'ATOM':
					for residue in residuesDict:
						if int(line[22:26]) == residue:
							residuesDict[residue].append(line)
							break
		# Create Residue object from the dictionary
		self.residues = {}
		for residue in residuesDict:
			self.residues[residue] = Residue(residuesDict[residue])