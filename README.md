# ProLIF
Protein-Ligand Interaction Fingerprints

## :small_blue_diamond: Description

ProLIF is a tool designed to generate Interaction FingerPrints (IFP) and compute similarity scores for protein-ligand interactions, given a reference ligand and a list of binding-site residues.

## :small_blue_diamond: Installation

ProLIF is written in Python 3, and uses the following non-standard libraries:
* [rdkit](http://www.rdkit.org/docs/Install.html)

Once all necessary libraries are installed, open a terminal and launch ProLIF with the following command:
```
python path/to/prolif.py --help
```
For Linux users, you can also make an alias for ease of use: add the following line in the `.bashrc` file located in your home directory: 
```
alias prolif="python path/to/prolif.py"
```

## :small_blue_diamond: Usage

```
python prolif.py --help
usage: prolif.py [-h] -r fileName -l fileName [fileName ...] -p fileName
                 --residues int [int ...] [--chain letter] [-o filename] [-v]
                 [--similarity {1,2,3,4,5,6,7,8,9,10,11,12}] [--alpha int]
                 [--beta int]

ProLIF: Protein Ligand Interaction Fingerprints
Generates Interaction FingerPrints (IFP) and a similarity score for protein-ligand interactions

optional arguments:
  -h, --help            show this help message and exit

INPUT arguments:
  -r fileName, --reference fileName
                        Path to your reference ligand. Type: SDF, MOL2, PDB
  -l fileName [fileName ...], --ligand fileName [fileName ...]
                        Path to your ligand(s). Type: SDF, MOL2, PDB
  -p fileName, --protein fileName
                        Path to your protein or bindingSite. Type: PDB
  --residues int [int ...]
                        Residues chosen for the interactions. Example: 11 32 34 171 174
  --chain letter        Restricts the IFP to the residues present in specified chain.

OUTPUT arguments:
  -o filename, --output filename
                        Path to the output CSV file
  -v, --verbose         Increase terminal output verbosity

Other arguments:
  --similarity {1,2,3,4,5,6,7,8,9,10,11,12}
                        Similarity score between molecule A and B :
                        Let 'a' and 'b' be the number of bits activated in molecules A and B,
                        'c' the number of activated bits in common, and 'd' the number of inactivated bits in common.
                            1) Tanimoto      : c/(a+b-c). Used by default
                            2) Dice          : 2c/(a+b)
                            3) Cosine        : c/sqrt(a*b)
                            4) Sokal         : c/(2a+2b-3c)
                            5) Russel        : c/(a+b+d-c)
                            6) RogotGoldberg : (c/(a+b))+(d/(a+b+2d-2c))
                            7) AllBit        : (a - (a XOR b))/a
                            8) Kulczynski    : c*(a+b)/(2*a*b)
                            9) McConnaughey  : (c*(a+b)-(a*b))/(a*b)
                           10) Asymmetric    : c/min(a,b)
                           11) BraunBlanquet : c/max(a,b)
                           12) Tversky       : c/(alpha*(a-c)+beta*(b-c)+c)
  --alpha int           Alpha parameter for Tversky
  --beta int            Beta parameter for Tversky

Mandatory arguments: --reference --ligand --protein --residues
1 single molecule per file.
```

## :small_blue_diamond: License

Unless otherwise noted, all files in this directory and all subdirectories are distributed under the Apache License, Version 2.0:
```
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
```

### TODO

- [ ] Put the IFP related methods in Ligand instead of Residue and Protein, and modify the main code accordingly
- [ ] Run with multithreading
- [ ] Evaluate the interactions (instead of returning randomly 0 or 1...)
