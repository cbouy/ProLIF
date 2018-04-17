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
from .residue import *
from .utils import *

class Protein:
    """Class for a protein"""

    def __init__(self, inputFile, reference, cutoff, residueList):
        """Initialization of the protein, defined by a list of residues"""
        self.residueList = residueList
        self.residues = {}
        self.inputFile = inputFile
        fileExtension = os.path.splitext(inputFile)[1]
        if fileExtension == '.mol2':
            self.residuesFromMOL2File()
        else:
            raise ValueError('{} files are not supported for the protein.'.format(fileExtension[1:].upper()))
        if self.residueList == None:
            self.residueList = self.detectCloseResidues(reference, cutoff)
        self.cleanResidues()

    def __repr__(self):
        return self.inputFile

    def residuesFromMOL2File(self):
        """Read a MOL2 file and assign each line to an object of class Atom"""
        # Create a molecule with RDKIT
        self.mol = Chem.MolFromMol2File(self.inputFile, sanitize=True, removeHs=False)
        # Read the atoms directly from the MOL2 file
        rec = mol2_reader(self.inputFile, ignoreH=False)
        # Loop through each RDKIT atom and assign them to a residue
        residues = {}
        coordinates = self.mol.GetConformer().GetPositions()
        for i,(rd_at, at) in enumerate(zip(self.mol.GetAtoms(), rec)):
            residue = at['residue']
            if residue not in residues:
                residues[residue] = []
            atom = {
                'resid'  :     residue,
                'resname':     residue[:3],
                'chain'  :     at['chain'],
                'coordinates': coordinates[i],
                'charge':      float(rd_at.GetProp('_TriposPartialCharge')),
                'aromatic':    rd_at.GetIsAromatic(),
                'ring':        rd_at.IsInRing(),
            }
            residues[residue].append(atom)
        # Create residues objects
        for residue in residues:
            self.residues[residue] = Residue(residues[residue])

    def detectCloseResidues(self, reference, cutoff):
        """Detect residues close to a reference ligand"""
        residueList = []
        for residue in self.residues:
            d = euclidianDistance(reference.centroid, self.residues[residue].centroid)
            if d <= cutoff:
                residueList.append(self.residues[residue].resid)
        return residueList

    def cleanResidues(self):
        """Cleans the residues of the protein to only keep those in self.residueList"""
        residues = {}
        for residue in self.residues:
            if residue in self.residueList:
                residues[residue] = self.residues[residue]
        self.residues = residues
