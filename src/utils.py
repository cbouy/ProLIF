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

from math import sqrt
import re

def euclidianDistance(a,b):
    """Euclidian distance between 2 lists of 3D coordinates"""
    return sqrt(sum([(xa-xb)**2 for xa,xb in zip(a,b)]))

def getCentroid(atoms):
    """Centroid coordinates (XYZ) of a list of atoms"""
    return [sum([atom['coordinates'][i] for atom in atoms])/len(atoms) for i in range(3)]

def mol2_reader(mol2_file, ignoreH):
    '''A simple MOL2 file reader. Can only read files with a single molecule in it.
    Returns a molecules. Each molecule is a list of atoms, where an atom
    is a dictionnary containing atom informations.'''
    num_atoms_lines = []
    first_lines     = []

    # Read file
    with open(mol2_file, "r") as f:
        lines = f.readlines()

    # Search for the line where the number of atoms is, and the first line where atom coordinates are readable
    for i, line in enumerate(lines):
        search_molecule = re.search(r'@<TRIPOS>MOLECULE', line)
        search_atom     = re.search(r'@<TRIPOS>ATOM', line)

        if search_molecule:
            # line with the number of atoms
            num_atoms_lines.append(i+2)
        elif search_atom:
            # first line with atom coordinates
            first_lines.append(i+1)

    for num_atoms_line, first_line in zip(num_atoms_lines, first_lines):
        molecule = get_mol_from_mol2(num_atoms_line, first_line, lines, ignoreH)

    return molecule

def get_mol_from_mol2(num_atoms_line, first_line, lines, ignoreH):
    '''Extracts a molecule from a mol2 file.
    num_atoms_line: index of the line containing the number of atoms, bonds...etc.
    first_line: index of the first line of the molecule to be extracted
    Returns a molecule as a list of atoms. An atom is a dictionnary containing
    atom informations that aren't read by the RDkit reader: index, chain, residue'''
    # Read number of atoms directly from the corresponding line
    data      = lines[num_atoms_line].split()
    num_atoms = int(data[0])

    # create molecule containing a list of atoms
    MOL = []
    # Fill the list with atomic parameters and coordinates
    for line in range(first_line, first_line + num_atoms):
        data = lines[line].split()
        if ignoreH: # if ignore H
            if data[5] == 'H': # if the atom read is an H
                continue # skip this atom
        # if it's not an H, or we don't ignore H, add this atom
        MOL.append(
        {
            'chain'   : int(data[6]),
            'residue' : str(data[7]),
        })
    return MOL
