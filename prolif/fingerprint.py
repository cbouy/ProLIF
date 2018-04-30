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
from rdkit import Chem, DataStructs
from .prolif import logger
from .utils import (getCentroid, euclidianDistance, getAngle, pointsToVector,
                    isinAngleLimits, getNormalVector)

class Fingerprint:
    """Class that generates an interaction fingerprint between a protein and a ligand"""

    def __init__(self,
        json_file=path.join(path.dirname(__file__),'parameters.json'),
        interactions=['HBdonor','HBacceptor','cation','anion','FaceToFace','EdgeToFace','hydrophobic']):
        # read parameters from json file
        with open(json_file) as data_file:
            logger.debug('Reading JSON parameters file from {}'.format(json_file))
            self.prm = json.load(data_file)
        # create aromatic patterns from json file
        self.AROMATIC_PATTERNS = [ Chem.MolFromSmarts(smart) for smart in self.prm["aromatic"]["smarts"]]
        # read interactions to compute
        self.interactions = interactions
        logger.info('Built fingerprint generator using the following bitstring: {}'.format(' '.join(self.interactions)))


    def __repr__(self):
        return ' '.join(self.interactions)


    def hasHydrophobic(self, ligand, residue):
        """Get the presence or absence of an hydrophobic interaction between
        a residue and a ligand"""
        hydrophobic = Chem.MolFromSmarts(self.prm["hydrophobic"]["smarts"])
        lig_matches = ligand.mol.GetSubstructMatches(hydrophobic)
        res_matches = residue.mol.GetSubstructMatches(hydrophobic)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    dist = euclidianDistance(
                        ligand.coordinates[lig_match[0]],
                        residue.coordinates[res_match[0]])
                    if dist <= self.prm["hydrophobic"]["distance"]:
                        return 1
        return 0


    def hasHBdonor(self, ligand, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as an acceptor and a ligand as a donor"""
        donor    = Chem.MolFromSmarts(self.prm["HBond"]["donor"])
        acceptor = Chem.MolFromSmarts(self.prm["HBond"]["acceptor"])
        lig_matches = ligand.mol.GetSubstructMatches(donor)
        res_matches = residue.mol.GetSubstructMatches(acceptor)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    # A-H ... B
                    a = ligand.coordinates[lig_match[0]]
                    h = ligand.coordinates[lig_match[1]]
                    b = residue.coordinates[res_match[0]]
                    dist = euclidianDistance(a, b)
                    if dist <= self.prm["HBond"]["distance"]:
                        ha = pointsToVector(h, a)
                        hb = pointsToVector(h, b)
                        angle = getAngle(ha, hb)
                        if isinAngleLimits(angle, *self.prm["HBond"]["angle"]):
                            return 1
        return 0


    def hasHBacceptor(self, ligand, residue):
        """Get the presence or absence of a H-bond interaction between
        a residue as a donor and a ligand as an acceptor"""
        donor    = Chem.MolFromSmarts(self.prm["HBond"]["donor"])
        acceptor = Chem.MolFromSmarts(self.prm["HBond"]["acceptor"])
        lig_matches = ligand.mol.GetSubstructMatches(acceptor)
        res_matches = residue.mol.GetSubstructMatches(donor)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    # A-H ... B
                    a = residue.coordinates[res_match[0]]
                    h = residue.coordinates[res_match[1]]
                    b = ligand.coordinates[lig_match[0]]
                    dist = euclidianDistance(a, b)
                    if dist <= self.prm["HBond"]["distance"]:
                        ha = pointsToVector(h, a)
                        hb = pointsToVector(h, b)
                        angle = getAngle(ha, hb)
                        if isinAngleLimits(angle, *self.prm["HBond"]["angle"]):
                            return 1
        return 0


    def hasXBdonor(self, ligand, residue):
        """Get the presence or absence of a Halogen Bond where the ligand acts as
        a donor"""
        donor    = Chem.MolFromSmarts(self.prm["XBond"]["donor"])
        acceptor = Chem.MolFromSmarts(self.prm["XBond"]["acceptor"])
        lig_matches = ligand.mol.GetSubstructMatches(donor)
        res_matches = residue.mol.GetSubstructMatches(acceptor)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    # A-X ... B
                    a = ligand.coordinates[lig_match[0]]
                    x = ligand.coordinates[lig_match[1]]
                    b = residue.coordinates[res_match[0]]
                    dist = euclidianDistance(a, b)
                    if dist <= self.prm["XBond"]["distance"]:
                        xa = pointsToVector(x, a)
                        xb = pointsToVector(x, b)
                        angle = getAngle(xa, xb)
                        if isinAngleLimits(angle, *self.prm["XBond"]["angle"]):
                            return 1
        return 0


    def hasXBacceptor(self, ligand, residue):
        """Get the presence or absence of a Halogen Bond where the residue acts as
        a donor"""
        donor    = Chem.MolFromSmarts(self.prm["XBond"]["donor"])
        acceptor = Chem.MolFromSmarts(self.prm["XBond"]["acceptor"])
        lig_matches = ligand.mol.GetSubstructMatches(acceptor)
        res_matches = residue.mol.GetSubstructMatches(donor)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    # A-X ... B
                    a = residue.coordinates[res_match[0]]
                    x = residue.coordinates[res_match[1]]
                    b = ligand.coordinates[lig_match[0]]
                    dist = euclidianDistance(a, b)
                    if dist <= self.prm["XBond"]["distance"]:
                        xa = pointsToVector(x, a)
                        xb = pointsToVector(x, b)
                        angle = getAngle(xa, xb)
                        if isinAngleLimits(angle, *self.prm["XBond"]["angle"]):
                            return 1
        return 0


    def hasCationic(self, ligand, residue):
        """Get the presence or absence of an ionic interaction between a residue
        as an anion and a ligand as a cation"""
        cation = Chem.MolFromSmarts(self.prm["ionic"]["cation"])
        anion  = Chem.MolFromSmarts(self.prm["ionic"]["anion"])
        lig_matches = ligand.mol.GetSubstructMatches(cation)
        res_matches = residue.mol.GetSubstructMatches(anion)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    dist = euclidianDistance(
                        ligand.coordinates[lig_match[0]],
                        residue.coordinates[res_match[0]])
                    if dist <= self.prm["ionic"]["distance"]:
                        return 1
        return 0


    def hasAnionic(self, ligand, residue):
        """Get the presence or absence of an ionic interaction between a residue
        as a cation and a ligand as an anion"""
        cation = Chem.MolFromSmarts(self.prm["ionic"]["cation"])
        anion  = Chem.MolFromSmarts(self.prm["ionic"]["anion"])
        lig_matches = ligand.mol.GetSubstructMatches(anion)
        res_matches = residue.mol.GetSubstructMatches(cation)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    dist = euclidianDistance(
                        ligand.coordinates[lig_match[0]],
                        residue.coordinates[res_match[0]])
                    if dist <= self.prm["ionic"]["distance"]:
                        return 1
        return 0


    def hasPiCation(self, ligand, residue):
        """Get the presence or absence of an interaction between a residue as
        a cation and a ligand as a pi system"""
        cation = Chem.MolFromSmarts(self.prm["ionic"]["cation"])
        res_matches = residue.mol.GetSubstructMatches(cation)
        if res_matches:
            for pi in self.AROMATIC_PATTERNS:
                lig_matches = ligand.mol.GetSubstructMatches(pi)
                for res_match in res_matches:
                    for lig_match in lig_matches:
                        pi_coords = ligand.coordinates[list(lig_match)]
                        centroid  = getCentroid(pi_coords)
                        resid_coords = residue.coordinates[res_match[0]]
                        dist = euclidianDistance(resid_coords, centroid)
                        if dist <= self.prm["pi-cation"]["distance"]:
                            # get vector in the ring plane
                            pi_vector = pointsToVector(centroid, pi_coords[0])
                            # vector normal to the ring plane
                            normal_pi = getNormalVector(pi_vector)
                            # vector between the centroid and the charge
                            centroid_charge = pointsToVector(centroid, resid_coords)
                            angle = getAngle(normal_pi, centroid_charge)
                            if isinAngleLimits(angle, *self.prm["pi-cation"]["angle"]):
                                return 1
        return 0


    def hasCationPi(self, ligand, residue):
        """Get the presence or absence of an interaction between a residue as a
        pi system and a ligand as a cation"""
        cation = Chem.MolFromSmarts(self.prm["ionic"]["cation"])
        lig_matches = ligand.mol.GetSubstructMatches(cation)
        if lig_matches:
            for pi in self.AROMATIC_PATTERNS:
                res_matches = residue.mol.GetSubstructMatches(pi)
                for res_match in res_matches:
                    for lig_match in lig_matches:
                        pi_coords = residue.coordinates[list(res_match)]
                        centroid  = getCentroid(pi_coords)
                        lig_coords = ligand.coordinates[lig_match[0]]
                        dist = euclidianDistance(lig_coords, centroid)
                        if dist <= self.prm["pi-cation"]["distance"]:
                            # vector between the centroid and a point of the plane
                            pi_vector = pointsToVector(centroid, pi_coords[0])
                            # normal vector
                            normal_pi = getNormalVector(pi_vector)
                            # vector between the centroid and the charge
                            centroid_charge = pointsToVector(centroid, lig_coords)
                            angle = getAngle(normal_pi, centroid_charge)
                            if isinAngleLimits(angle, *self.prm["pi-cation"]["angle"]):
                                return 1
        return 0


    def hasFaceToFace(self, ligand, residue):
        """Get the presence or absence of an aromatic face to face interaction
        between a residue and a ligand"""
        for pi_res in self.AROMATIC_PATTERNS:
            for pi_lig in self.AROMATIC_PATTERNS:
                res_matches = residue.mol.GetSubstructMatches(pi_res)
                lig_matches = ligand.mol.GetSubstructMatches(pi_lig)
                for lig_match in lig_matches:
                    for res_match in res_matches:
                        res_pi_coords = residue.coordinates[list(res_match)]
                        res_centroid  = getCentroid(res_pi_coords)
                        lig_pi_coords = ligand.coordinates[list(lig_match)]
                        lig_centroid  = getCentroid(lig_pi_coords)
                        dist = euclidianDistance(lig_centroid, res_centroid)
                        if dist <= self.prm["aromatic"]["distance"]:
                            # ligand
                            pi_vector = pointsToVector(lig_centroid, lig_pi_coords[0])
                            lig_normal = getNormalVector(pi_vector)
                            # residue
                            pi_vector = pointsToVector(res_centroid, res_pi_coords[0])
                            res_normal = getNormalVector(pi_vector)
                            # angle
                            angle = getAngle(res_normal, lig_normal)
                            if isinAngleLimits(angle, *self.prm["aromatic"]["FaceToFace"]):
                                return 1
        return 0


    def hasEdgeToFace(self, ligand, residue):
        """Get the presence or absence of an aromatic face to edge interaction
        between a residue and a ligand"""
        for pi_res in self.AROMATIC_PATTERNS:
            for pi_lig in self.AROMATIC_PATTERNS:
                res_matches = residue.mol.GetSubstructMatches(pi_res)
                lig_matches = ligand.mol.GetSubstructMatches(pi_lig)
                for lig_match in lig_matches:
                    for res_match in res_matches:
                        res_pi_coords = residue.coordinates[list(res_match)]
                        res_centroid  = getCentroid(res_pi_coords)
                        lig_pi_coords = ligand.coordinates[list(lig_match)]
                        lig_centroid  = getCentroid(lig_pi_coords)
                        dist = euclidianDistance(lig_centroid, res_centroid)
                        if dist <= self.prm["aromatic"]["distance"]:
                            # ligand
                            pi_vector = pointsToVector(lig_centroid, lig_pi_coords[0])
                            lig_normal = getNormalVector(pi_vector)
                            # residue
                            pi_vector = pointsToVector(res_centroid, res_pi_coords[0])
                            res_normal = getNormalVector(pi_vector)
                            # angle
                            angle = getAngle(res_normal, lig_normal)
                            if isinAngleLimits(angle, *self.prm["aromatic"]["EdgeToFace"]):
                                return 1
        return 0


    def hasMetalDonor(self, ligand, residue):
        """Get the presence or absence of a metal complexation where the ligand is a metal"""
        metal = Chem.MolFromSmarts(self.prm["metallic"]["metal"])
        lig   = Chem.MolFromSmarts(self.prm["metallic"]["ligand"])
        lig_matches = ligand.mol.GetSubstructMatches(metal)
        res_matches = residue.mol.GetSubstructMatches(lig)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    dist = euclidianDistance(
                        ligand.coordinates[lig_match[0]],
                        residue.coordinates[res_match[0]])
                    if dist <= self.prm["metallic"]["distance"]:
                        return 1
        return 0


    def hasMetalAcceptor(self, ligand, residue):
        """Get the presence or absence of a metal complexation where the residue is a metal"""
        metal = Chem.MolFromSmarts(self.prm["metallic"]["metal"])
        lig   = Chem.MolFromSmarts(self.prm["metallic"]["ligand"])
        lig_matches = ligand.mol.GetSubstructMatches(lig)
        res_matches = residue.mol.GetSubstructMatches(metal)
        if lig_matches and res_matches:
            for lig_match in lig_matches:
                for res_match in res_matches:
                    dist = euclidianDistance(
                        ligand.coordinates[lig_match[0]],
                        residue.coordinates[res_match[0]])
                    if dist <= self.prm["metallic"]["distance"]:
                        return 1
        return 0


    def generateBitstring(self, ligand, residue):
        """Generate the complete bitstring for the interactions of a residue with a ligand"""
        bitstring = []
        for interaction in self.interactions:
            if   interaction == 'HBdonor':
                bitstring.append(self.hasHBdonor(ligand, residue))
            elif interaction == 'HBacceptor':
                bitstring.append(self.hasHBacceptor(ligand, residue))
            elif interaction == 'XBdonor':
                bitstring.append(self.hasXBdonor(ligand, residue))
            elif interaction == 'XBacceptor':
                bitstring.append(self.hasXBdonor(ligand, residue))
            elif interaction == 'cation':
                bitstring.append(self.hasCationic(ligand, residue))
            elif interaction == 'anion':
                bitstring.append(self.hasAnionic(ligand, residue))
            elif interaction == 'FaceToFace':
                bitstring.append(self.hasFaceToFace(ligand, residue))
            elif interaction == 'EdgeToFace':
                bitstring.append(self.hasEdgeToFace(ligand, residue))
            elif interaction == 'pi-cation':
                bitstring.append(self.hasPiCation(ligand, residue))
            elif interaction == 'cation-pi':
                bitstring.append(self.hasCationPi(ligand, residue))
            elif interaction == 'hydrophobic':
                bitstring.append(self.hasHydrophobic(ligand, residue))
            elif interaction == 'MBdonor':
                bitstring.append(self.hasMetalDonor(ligand, residue))
            elif interaction == 'MBacceptor':
                bitstring.append(self.hasMetalAcceptor(ligand, residue))
        return ''.join(str(bit) for bit in bitstring)


    def generateIFP(self, ligand, protein):
        """Generates the complete IFP from each residue's bitstring"""
        IFPvector = DataStructs.ExplicitBitVect(len(self.interactions)*len(protein.residues))
        i=0
        IFP = ''
        for residue in sorted(protein.residues, key=lambda x: x[3:]):
            bitstring = self.generateBitstring(ligand, protein.residues[residue])
            for bit in bitstring:
                if bit == '1':
                    IFPvector.SetBit(i)
                i+=1
            IFP += bitstring
        ligand.setIFP(IFP, IFPvector)
