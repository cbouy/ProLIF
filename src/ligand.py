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

from rdkit import Chem, DataStructs
import os.path

class Ligand:
    """Class for a ligand"""
    def __init__(self, inputFile):
        """Initialize the ligand from a file"""
        self.inputFile = inputFile
        fileExtension = os.path.splitext(inputFile)[1]
        if fileExtension == '.mol2':
            self.mol = Chem.MolFromMol2File(inputFile, sanitize=True, removeHs=False)
        else:
            raise ValueError('{} files are not supported for the ligand.'.format(fileExtension[1:].upper()))
        # Set Centroid
        coordinates = self.mol.GetConformer().GetPositions()
        self.centroid = [sum([atom[i] for atom in coordinates])/len(coordinates) for i in range(3)]

    def hasHydrophobic(self, residue):
        """Get the presence or absence of an hydrophobic interaction between
        a residue and a ligand"""
        return 0

    def hasArFtoF(self, residue):
        """Get the presence or absence of an aromatic face to face interaction
        between a residue and a ligand"""
        return 0

    def hasArFtoE(self, residue):
        """Get the presence or absence of an aromatic face to edge interaction
        between a residue and a ligand"""
        return 0

    def hasHBd(self, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as an acceptor and a ligand as a donor"""
        return 0

    def hasHBa(self, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as a donor and a ligand as an acceptor"""
        return 0

    def hasXB(self, residue):
        """Get the presence or absence of a Halogen Bond"""
        return 0

    def hasCationic(self, residue):
        """Get the presence or absence of an ionic interaction between a residue and a ligand as a cation"""
        return 0

    def hasAnionic(self, residue):
        """Get the presence or absence of an ionic interaction between a residue and a ligand as an anion"""
        return 0

    def hasPiCation(self, residue):
        """Get the presence or absence of Pi-Cation interaction"""
        return 0

    def hasMetalComplexation(self, residue):
        """Get the presence or absence of a metal complexation"""
        return 0

    def generateBitstring(self, residue):
        """Generate the complete bitstring for the interactions of a residue with a ligand"""
        estimate = [
            self.hasHydrophobic(residue),
            self.hasArFtoF(residue),
            self.hasArFtoE(residue),
            self.hasHBd(residue),
            self.hasHBa(residue),
            self.hasXB(residue),
            self.hasCationic(residue),
            self.hasAnionic(residue),
            ]
        return ''.join(str(bit) for bit in estimate)

    def generateIFP(self, protein):
        """Generates the complete IFP from each residue's bitstring"""
        self.IFPvector = DataStructs.ExplicitBitVect(8*len(protein.residues))
        i=0
        self.IFP = ''
        for residue in protein.residues:
            bitstring = self.generateBitstring(protein.residues[residue])
            for value in bitstring:
                if value == '1':
                    self.IFPvector.SetBit(i)
                i+=1
            self.IFP += bitstring

    def setSimilarity(self, S):
        """Set the value for the similarity score between the ligand and a reference"""
        self.S = S
