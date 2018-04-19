[![PyPI - Version](https://badge.fury.io/py/prolif.svg)](https://pypi.org/project/prolif/)
[![PyPI - License](https://img.shields.io/pypi/l/prolif.svg)](https://pypi.org/project/prolif/)
[![PyPI - Status](https://img.shields.io/pypi/status/prolif.svg)](https://pypi.org/project/prolif/)
[![Build Status](https://travis-ci.org/cbouy/ProLIF.svg?branch=master)](https://travis-ci.org/cbouy/ProLIF)
[![Coverage Status](https://coveralls.io/repos/github/cbouy/ProLIF/badge.svg?branch=master)](https://coveralls.io/github/cbouy/ProLIF?branch=master)

# ProLIF
Protein-Ligand Interaction Fingerprints

:warning: This project is under development, do not use it in the current state :warning:

## Description

ProLIF is a tool designed to generate Interaction FingerPrints (IFP) and compute similarity scores for protein-ligand interactions, given a reference ligand and a list of binding-site residues.

## Installation

ProLIF is written in Python 3, and uses the following non-standard libraries:
* [rdkit](http://www.rdkit.org/docs/Install.html)

To install rdkit with Anaconda, use the following command:
```
conda install -c rdkit rdkit
```

## Usage

```
INPUT arguments:
  -r fileName, --reference fileName
                        Path to your reference ligand.
  -l fileName [fileName ...], --ligand fileName [fileName ...]
                        Path to your ligand(s).
  -p fileName, --protein fileName
                        Path to your protein.
  --residues RESIDUES [RESIDUES ...]
                        Residues chosen for the interactions. Default: automatically detect residues within --cutoff of the reference ligand
  --cutoff float        Cutoff for automatic residue detection. Default: 5.0 angströms
  --json fileName       Path to a custom parameters file. Default: prolif.json

OUTPUT arguments:
  -o filename, --output filename
                        Path to the output CSV file
  -v, --verbose         Increase terminal output verbosity

Other arguments:
  --interactions bit [bit ...]
                        List of interactions used to build the fingerprint.
                        -) hydrogen bond: HBdonor, HBacceptor
                        -) halogen bond:  XBdonor
                        -) ionic: cation, anion
                        -) pi-stacking: FaceToFace, FaceToEdge
                        -) hydrophobic
                        -) pi-cation
                        -) metal
                        Default: HBdonor HBacceptor cation anion FaceToFace FaceToEdge pi-cation hydrophobic
  --score {tanimoto,dice,tversky}
                        Similarity score between molecule A and B :
                        Let 'a' and 'b' be the number of bits activated in molecules A and B, and 'c' the number of activated bits in common.
                        -) tanimoto : c/(a+b-c). Used by default
                        -) dice     : 2c/(a+b)
                        -) tversky  : c/(alpha*(a-c)+beta*(b-c)+c)
  --alpha int           Alpha parameter for Tversky. Default: 0.7
  --beta int            Beta parameter for Tversky. Default: 0.3

Mandatory arguments: --reference --ligand --protein
MOL2 files only.
```

## License

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
