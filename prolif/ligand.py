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

    def __repr__(self):
        return self.inputFile

    def setIFP(self, IFP, vector):
        """Set the IFP for the ligand, as a bitstring and vector"""
        self.IFP = IFP
        self.IFPvector = vector

    def getSimilarity(self, reference, method='tanimoto', alpha=None, beta=None):
        if   method == 'tanimoto':
            return DataStructs.TanimotoSimilarity(reference.IFPvector, self.IFPvector)
        elif method == 'dice':
            return DataStructs.DiceSimilarity(reference.IFPvector, self.IFPvector)
        elif method == 'tversky':
            return DataStructs.TverskySimilarity(reference.IFPvector, self.IFPvector, alpha, beta)

    def setSimilarity(self, score):
        """Set the value for the similarity score between the ligand and a reference"""
        self.score = score
