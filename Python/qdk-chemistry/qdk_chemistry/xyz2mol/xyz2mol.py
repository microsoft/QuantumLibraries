"""
Module for generating rdkit molobj/smiles/molecular graph from free atoms
Implementation by Jan H. Jensen, based on the paper
    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
    to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
    DOI: 10.1002/bkcs.10334
"""

import copy
import itertools
import sys

from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
try:
    from rdkit.Chem import rdEHTTools #requires RDKit 2019.9.1 or later
except ImportError:
    rdEHTTools = None

from collections import defaultdict

import numpy as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

from qdk_chemistry.xyz2mol.ac import xyz2AC, AC2mol
from qdk_chemistry.xyz2mol.util import chiral_stereo_check


def xyz2mol(
    atoms, 
    coordinates,
    charge=0,
    allow_charged_fragments=True,
    use_graph=True,
    use_huckel=False,
    embed_chiral=True
):
    """
    Generate a rdkit molobj from atoms, coordinates and a total_charge.
    args:
        atoms - list of atom types (int)
        coordinates - 3xN Cartesian coordinates
        charge - total charge of the system (default: 0)
    optional:
        allow_charged_fragments - alternatively radicals are made
        use_graph - use graph (networkx)
        use_huckel - Use Huckel method for atom connectivity prediction
        embed_chiral - embed chiral information to the molecule
    returns:
        mols - list of rdkit molobjects
    """

    # Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge,
    # and mol object with no connectivity information
    AC, mol = xyz2AC(atoms, coordinates, charge, use_huckel=use_huckel)

    # Convert AC to bond order matrix and add connectivity and charge info to
    # mol object
    new_mols = AC2mol(mol, AC, atoms, charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph)

    # Check for stereocenters and chiral centers
    if embed_chiral:
        for new_mol in new_mols:
            chiral_stereo_check(new_mol)

    return new_mols
