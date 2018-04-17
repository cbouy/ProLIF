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
from .utils import getCentroid

class Residue:
    """Class for a residue in a protein"""

    def __init__(self, atoms):
        self.atoms    = atoms
        self.resname  = atoms[0]['resname'] # 3 letter code
        self.resid    = atoms[0]['resid']   # unique identifier for the residue
        self.centroid = getCentroid(self.atoms)

    def __repr__(self):
        return self.resid
