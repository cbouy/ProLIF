#!/usr/bin/python3

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

from .ligand import *
from .protein import *
from .fingerprint import *

def main(args):

    # Read files
    fingerprint = Fingerprint(args.json, args.interactions)
    if args.verbose:
        print('Built fingerprint generator using the following bitstring:')
        print(fingerprint)
    reference   = Ligand(args.reference)
    if args.verbose:
        print('Read reference molecule from ' + str(reference))
    protein     = Protein(args.protein, reference, cutoff=args.cutoff, residueList=args.residues)
    if args.verbose:
        print('Read protein from ' + str(protein))
    length      = len(args.interactions)
    # Print residues on terminal:
    for residue in protein.residues:
        print('{resid: >{length}s}'.format(resid=protein.residues[residue].resid, length=length), end='')
    print()
    # Generate the IFP between the reference ligand and the protein
    fingerprint.generateIFP(reference, protein)
    print(reference.IFP)
    # Loop over ligands:
    ligandList = []
    for lig in args.ligand:
        ligand = Ligand(lig)
        # Generate the IFP between a ligand and the protein
        fingerprint.generateIFP(ligand, protein)
        # Calculate similarity
        score = ligand.getSimilarity(reference, args.score, args.alpha, args.beta)
        ligand.setSimilarity(score)
        print(ligand.IFP, '{:.4f}'.format(ligand.score))
        ligandList.append(ligand)
    # Output
    if args.output:
        with open(args.output, 'w') as file:
            file.write('File,SimilarityScore')
            for residue in protein.residues:
                file.write(',{}'.format(protein.residues[residue].resid))
            file.write('\n')
            CSIFP = ','.join(reference.IFP[i:i+length] for i in range(0, len(reference.IFP), length))
            file.write('{},{:.4f},{}\n'.format(args.reference, 1, CSIFP))
            for ligand in ligandList:
                CSIFP = ','.join(ligand.IFP[i:i+length] for i in range(0, len(ligand.IFP), length))
                file.write('{},{:.4f},{}\n'.format(ligand, ligand.score, CSIFP))

if __name__ == '__main__':
    main(args)
