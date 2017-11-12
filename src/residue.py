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
		self.setResName('{}{}'.format(pdbInfo.GetResidueName(),pdbInfo.GetResidueNumber()))
	
	def setResName(self, residueName):
		"""Set the name of the residue"""
		self.resName = residueName