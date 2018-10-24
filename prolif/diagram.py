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

from rdkit.Chem import Draw, rdDepictor
import rdkit.Geometry as rdGeometry
from copy import deepcopy
import numpy as np


STYLE = {
    'HBdonor': {
        'color': '#FF4500', # orangered
        'dash': '8,8',
        'width': '2',
    },
    'HBacceptor': {
        'color': '#00BFFF', # deepskyblue
        'dash' : '8,8',
        'width': '2',
    },
    'cation': {
        'color': '#FF0000', # red
        'dash': '2,10',
        'width': '6',
    },
    'anion': {
        'color': '#4169E1', # royalblue
        'dash': '2,10',
        'width': '6',
    },
    'XBdonor': {
        'color': '#FF00FF', # magenta
        'dash': '6,6',
        'width': '2',
    },
    'XBacceptor': {
        'color': '#9400D3', # darkviolet
        'dash': '6,6',
        'width': '2',
    },
    'hydrophobic': {
        'color': '#00FF00', # lime
        'dash': '1,5',
        'width': '8',
    },
    'FaceToFace':{
        'color': '#8B4513', # saddlebrown
        'dash': '6,4,2,4',
        'width': '2',
    },
    'EdgeToFace':{
        'color': '#7FFF00', # chartreuse
        'dash': '6,4,2,4',
        'width': '2',
    },
    'pi-cation':{
        'color': '#FFD700', # gold
        'dash': '6,4,2,4',
        'width': '2',
    },
    'cation-pi':{
        'color': '#00FFFF', # cyan
        'dash': '6,4,2,4',
        'width': '2',
    },
    'MBdonor':{
        'color': '#DC143C', # crimson
        'dash': '2,5',
        'width': '3',
    },
    'MBacceptor':{
        'color': '#6A5ACD', # slateblue
        'dash': '2,5',
        'width': '3',
    },
}

def generate_interaction_diagram(ligand, fingerprint, output=None, init=False):
    """Generate a SVG image depicting the interactions of the ligand with the
    residues of the protein"""
    # Create canvas with options
    drawer = Draw.MolDraw2DSVG(400,400)
    options = drawer.drawOptions()
    #options.useBWAtomPalette()
    options.padding = 0.2
    options.additionalAtomLabelPadding = 0.1
    # Draw molecule
    lig = deepcopy(ligand.mol)
    rdDepictor.GenerateDepictionMatching3DStructure(lig, ligand.mol)
    drawer.DrawMolecule(lig)
    drawer.FinishDrawing()
    # Separate the header of the SVG file from the molecule drawing
    img = drawer.GetDrawingText().replace('svg:','').replace('</svg>','')
    img = img.split('\n')
    header = img[:8]
    mol = img[8:]
    # Draw interactions
    interactions_drawing = []
    for resname in fingerprint.diagram[ligand.inputFile]:
        draw = []
        # TODO: replace by user drag-and-drop positions for the reference (use init=True)
        res_x = np.random.randint(20,380)
        res_y = np.random.randint(20,380)
        for interaction in fingerprint.diagram[ligand.inputFile][resname]:
            style = STYLE[interaction]
            # coordinates of the atom on the canvas
            if interaction == 'pi-cation': # coordinates of the centroid of the ring
                pi_atoms = fingerprint.diagram[ligand.inputFile][resname][interaction]
                x,y = 0,0
                for atom in pi_atoms:
                    at = drawer.GetDrawCoords(atom)
                    x += at.x
                    y += at.y
                at = rdGeometry.Point2D(x/len(pi_atoms), y/len(pi_atoms))
            else: # coordinates of the atom
                ligand_atom = fingerprint.diagram[ligand.inputFile][resname][interaction]
                at = drawer.GetDrawCoords(ligand_atom)
            # draw dashed line between interacting atom and residue
            draw += ['''<path d='M {rx},{ry} {ax},{ay}' stroke-dasharray='{dash}' style='fill:none;fill-rule:evenodd;stroke:{c};stroke-width:{w}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:0.7' />'''.format(
                c=style['color'], dash=style['dash'], w=style['width'],
                rx=res_x, ry=res_y, ax=at.x, ay=at.y)]
        # draw the residue
        draw += ['''\
<ellipse cx='{cx}' cy='{cy}' rx='{rx}' ry='12' style='fill:white;fill-rule:evenodd;stroke:black;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='{tx}' y='{ty}' style='font-size:13px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' ><tspan>{text}</tspan></text>'''.format(
            text=resname, tx=res_x-4*len(resname), ty=res_y+4,
            cx=res_x, cy=res_y, rx=int(len(resname)*5.8)-1)]
        # Add interactions to the
        interactions_drawing += draw
    img = '\n'.join(header + interactions_drawing + mol) + '</svg>'

    # Write file
    if not output:
        output = '{}.svg'.format(ligand.inputName)
    with open(output,'w') as f:
        f.write(img)
