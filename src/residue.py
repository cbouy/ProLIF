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

import random
from rdkit import Chem

class Residue:
	"""Class for a residue in a protein"""
	
	def __init__(self, atomList):
		"""Initialization of the residue, defined by it's residue number"""
		# Make a rdkit molecule object for the residue 
		block = ''.join(line for line in atomList)
		self.structure = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=True)
		# Retrieve information from PDB file and use it to set the residue name for the output
		pdbInfo = self.structure.GetAtoms()[0].GetPDBResidueInfo()
		self.setResName('{}{}'.format(pdbInfo.GetResidueName(),pdbInfo.GetResidueNumber()))
	
	def setResName(self, residueName):
		"""Set the name of the residue"""
		self.resName = residueName
	
	def hasHydrophobic(self, ligand):
		"""Get the presence or absence of an hydrophobic interaction between a residue and a ligand"""
		return random.randint(0,1)
	
	def hasArFtoF(self, ligand):
		"""Get the presence or absence of an aromatic face to face interaction between a residue and a ligand"""
		return random.randint(0,1)
	
	def hasArFtoE(self, ligand):
		"""Get the presence or absence of an aromatic face to edge interaction between a residue and a ligand"""
		return random.randint(0,1)
	
	def hasHBd(self, ligand):
		"""Get the presence or absence of a H-bond interaction between a residue as a donor and a ligand"""
		return random.randint(0,1)
	
	def hasHBa(self, ligand):
		"""Get the presence or absence of a H-bond interaction between a residue as an acceptor and a ligand"""
		return random.randint(0,1)

	def hasXB(self, ligand):
		"""Get the presence or absence of a Halogen Bond"""
		return random.randint(0,1)
	
	def hasCationic(self, ligand):
		"""Get the presence or absence of an ionic interaction between a residue as a cation and a ligand"""
		return random.randint(0,1)
	
	def hasAnionic(self, ligand):
		"""Get the presence or absence of an ionic interaction between a residue as an anion and a ligand"""
		return random.randint(0,1)

	def generateBitstring(self, ligand):
		"""Generate the complete bitstring for the interactions of a residue with a ligand"""
		estimate = [ self.hasHydrophobic(ligand),
			self.hasArFtoF(ligand),
			self.hasArFtoE(ligand),
			self.hasHBd(ligand),
			self.hasHBa(ligand),
			self.hasXB(ligand),
			self.hasCationic(ligand),
			self.hasAnionic(ligand) ]
		self.bitstring = ''.join(str(bit) for bit in estimate)