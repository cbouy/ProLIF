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

class Residue:
	"""Class for a residue in a protein"""
	
	def __init__(self, atomList):
		"""Initialization of the residue, defined by it's residue number"""
		# Make a rdkit molecule object for the residue 
		block = ''.join(line for line in atomList)
		self.mol = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=True)
		# Retrieve information from PDB file and use it to set the residue name for the output
		pdbInfo = self.mol.GetAtoms()[0].GetPDBResidueInfo()
		self.aminoAcid = pdbInfo.GetResidueName()
		self.setResName('{}{}'.format(self.aminoAcid,pdbInfo.GetResidueNumber()))
		self.setProfile()
	
	def setResName(self, residueName):
		"""Set the name of the residue"""
		self.resName = residueName

	def setProfile(self):
		"""Set the profile of a residue based on the amino-acid"""
		# Hydrophobic
		if self.aminoAcid in ['ILE', 'VAL', 'LEU', 'PHE', 'MET', 'ALA', 'GLY', 'PRO', 'TRP', 'TYR']:
			self.isHydrophobic = True
		else:
			self.isHydrophobic = False
		# Aromatic
		if self.aminoAcid in ['PHE', 'TRP', 'TYR', 'HIS']:
			self.isAromatic = True
		else:
			self.isAromatic = False
		# HB and XB acceptor
		if self.aminoAcid in ['TYR', 'SER', 'THR', 'ASN', 'GLN', 'CYS']:
			self.isHBa = True
		else:
			self.isHBa = False
		# HB donor
		if self.aminoAcid in ['TYR', 'SER', 'THR', 'ASN', 'GLN', 'TRP', 'ARG']:
			self.isHBd = True
		else:
			self.isHBd = False
		# Anionic
		if self.aminoAcid in ['ASP', 'GLU']:
			self.isAnionic = True
		else:
			self.isAnionic = False
		# Cationic
		if self.aminoAcid in ['ARG', 'HIS', 'LYS']:
			self.isCationic = True
		else:
			self.isCationic = False