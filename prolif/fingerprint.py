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

import json
from os import path
from rdkit import DataStructs
from random import randint

class Fingerprint:
    """Class that generates an interaction fingerprint between a protein and a ligand"""

    def __init__(self,
        json_file=path.join(path.dirname(__file__),'prolif.json'),
        interactions=['HBdonor','HBacceptor','cation','anion','FaceToFace','FaceToEdge','pi-cation','hydrophobic']):
        # read parameters from json file
        with open(json_file) as data_file:
            self.prm = json.load(data_file)
        # read interactions to compute
        self.interactions = interactions

    def __repr__(self):
        return ' '.join(self.interactions)

    def hasHydrophobic(self, ligand, residue):
        """Get the presence or absence of an hydrophobic interaction between
        a residue and a ligand"""
        return randint(0, 1)

    def hasFaceToFace(self, ligand, residue):
        """Get the presence or absence of an aromatic face to face interaction
        between a residue and a ligand"""
        return randint(0, 1)

    def hasFaceToEdge(self, ligand, residue):
        """Get the presence or absence of an aromatic face to edge interaction
        between a residue and a ligand"""
        return randint(0, 1)

    def hasHBdonor(self, ligand, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as an acceptor and a ligand as a donor"""
        return randint(0, 1)

    def hasHBacceptor(self, ligand, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as a donor and a ligand as an acceptor"""
        return randint(0, 1)

    def hasXBdonor(self, ligand, residue):
        """Get the presence or absence of a Halogen Bond"""
        return randint(0, 1)

    def hasCationic(self, ligand, residue):
        """Get the presence or absence of an ionic interaction between a residue and a ligand as a cation"""
        return randint(0, 1)

    def hasAnionic(self, ligand, residue):
        """Get the presence or absence of an ionic interaction between a residue and a ligand as an anion"""
        return randint(0, 1)

    def hasPiCation(self, ligand, residue):
        """Get the presence or absence of Pi-Cation interaction"""
        return randint(0, 1)

    def hasMetal(self, ligand, residue):
        """Get the presence or absence of a metal complexation"""
        return randint(0, 1)

    def generateBitstring(self, ligand, residue):
        """Generate the complete bitstring for the interactions of a residue with a ligand"""
        bitstring = []
        if 'HBdonor' in self.interactions:
            bitstring.append(self.hasHBdonor(ligand, residue))
        if 'HBacceptor' in self.interactions:
            bitstring.append(self.hasHBacceptor(ligand, residue))
        if 'cation' in self.interactions:
            bitstring.append(self.hasCationic(ligand, residue))
        if 'anion' in self.interactions:
            bitstring.append(self.hasAnionic(ligand, residue))
        if 'FaceToFace' in self.interactions:
            bitstring.append(self.hasFaceToFace(ligand, residue))
        if 'FaceToEdge' in self.interactions:
            bitstring.append(self.hasFaceToEdge(ligand, residue))
        if 'pi-cation' in self.interactions:
            bitstring.append(self.hasPiCation(ligand, residue))
        if 'hydrophobic' in self.interactions:
            bitstring.append(self.hasHydrophobic(ligand, residue))
        if 'XBdonor' in self.interactions:
            bitstring.append(self.hasXBdonor(ligand, residue))
        if 'metal' in self.interactions:
            bitstring.append(self.hasMetal(ligand, residue))
        return ''.join(str(bit) for bit in bitstring)

    def generateIFP(self, ligand, protein):
        """Generates the complete IFP from each residue's bitstring"""
        IFPvector = DataStructs.ExplicitBitVect(len(self.interactions)*len(protein.residues))
        i=0
        IFP = ''
        for residue in protein.residues:
            bitstring = self.generateBitstring(ligand, protein.residues[residue])
            for value in bitstring:
                if value == '1':
                    IFPvector.SetBit(i)
                i+=1
            IFP += bitstring
        ligand.setIFP(IFP, IFPvector)
