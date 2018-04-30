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

import logging
from rdkit import RDLogger
from .logger import logger, stream_handler
from .ligand import *
from .protein import *
from .fingerprint import *

def main(args):
    # set logger level
    lg = RDLogger.logger()
    if args.log == 'CRITICAL':
        lg.setLevel(RDLogger.CRITICAL)
        stream_handler.setLevel(logging.CRITICAL)
    elif args.log == 'ERROR':
        lg.setLevel(RDLogger.ERROR)
        stream_handler.setLevel(logging.ERROR)
    elif args.log == 'WARNING':
        lg.setLevel(RDLogger.WARNING)
        stream_handler.setLevel(logging.WARNING)
    elif args.log == 'INFO':
        lg.setLevel(RDLogger.INFO)
        stream_handler.setLevel(logging.INFO)
    elif args.log == 'DEBUG':
        lg.setLevel(RDLogger.DEBUG)
        stream_handler.setLevel(logging.DEBUG)
    logger.info('Using {} to compute similarity between fingerprints'.format(args.score))
    # Read files
    fingerprint = Fingerprint(args.json, args.interactions)
    reference   = Ligand(args.reference)
    protein     = Protein(args.protein, reference, cutoff=args.cutoff, residueList=args.residues)
    length      = len(args.interactions)
    residues    = [protein.residues[residue] for residue in sorted(protein.residues, key=lambda x: int(x[3:]))]
    # Print residues on terminal:
    print(''.join('{resname: <{length}s}'.format(
        resname=residue.resname, length=length
        ) for residue in residues))
    # Generate the IFP between the reference ligand and the protein
    fingerprint.generateIFP(reference, protein)
    ifp_list = [reference.IFP[i:i+length] for i in range(0, len(reference.IFP), length)]
    print(''.join('{ifp: <{size}s}'.format(
        ifp=ifp_list[i], size=len(residues[i].resname)
        ) for i in range(len(ifp_list))), reference.inputFile)

    # Loop over ligands:
    ligandList = []
    for lig in args.ligand:
        ligand = Ligand(lig)
        # Generate the IFP between a ligand and the protein
        fingerprint.generateIFP(ligand, protein)
        # Calculate similarity
        score = ligand.getSimilarity(reference, args.score, args.alpha, args.beta)
        ligand.setSimilarity(score)
        ifp_list = [ligand.IFP[i:i+length] for i in range(0, len(ligand.IFP), length)]
        print(''.join('{ifp: <{size}s}'.format(ifp=ifp_list[i], size=len(residues[i].resname)) for i in range(len(ifp_list)) ),
              '{:.3f}'.format(ligand.score), ligand.inputFile)
        ligandList.append(ligand)

    # Output
    if args.output:
        logger.info('Writing CSV formatted output to ' + args.output)
        with open(args.output, 'w') as f:
            f.write('File,SimilarityScore')
            for residue in residues:
                f.write(',{}'.format(protein.residues[residue].resname))
            f.write('\n')
            CSIFP = ','.join(reference.IFP[i:i+length] for i in range(0, len(reference.IFP), length))
            f.write('{},,{}\n'.format(args.reference, CSIFP))
            for ligand in ligandList:
                CSIFP = ','.join(ligand.IFP[i:i+length] for i in range(0, len(ligand.IFP), length))
                f.write('{},{:.3f},{}\n'.format(ligand.inputFile, ligand.score, CSIFP))

if __name__ == '__main__':
    main(args)
