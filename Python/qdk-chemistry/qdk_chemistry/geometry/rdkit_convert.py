# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for converting to and from RDKit data structures
"""

from typing import List, Tuple, Iterable, TYPE_CHECKING

from .xyz import coordinates_to_xyz

from rdkit.Chem import AllChem as Chem

if TYPE_CHECKING:
    from rdkit.Chem import Mol, Conformer


def get_conformer(mol: "Mol", num_confs: int=10) -> "Conformer":
    """Get lowest-energy Conformer for molecular fragment.
    If conformers don't converge, get lowest energy conformer.

    :param mol: Molecular fragment
    :type mol: Mol
    :param num_confs: Number of configurations to generate, defaults to 10.
    :type num_confs: int

    :return: Lowest-energy conformer for molecule
    :rtype: Conformer
    """
    conformers = mol.GetConformers()
    if len(conformers) < num_confs:
        # Embed num_confs conformers into molecule object
        Chem.EmbedMultipleConfs(mol, numConfs=num_confs)
        conformers = mol.GetConformers()

    res = Chem.MMFFOptimizeMoleculeConfs(mol, numThreads=0) # resturns a list of tuples (not_converged, energy)
    if 0 in [not_converged for not_converged, _ in res]:
        idx = res.index(min(res, key=lambda x: x[1]))
        conformer = conformers[idx]
    else:
        ((not_converged, energy), conformer) = min(zip(res, mol.GetConformers()), key=lambda x: x[0][1])
        if not_converged:
            print(f"Solution did not converge. Lowest energy found: {energy}")
    return conformer


def _conformer_to_coordinates(
        symbols: Iterable[str], 
        conformer: "Conformer"
    ) -> List[Tuple[str, float, float, float]]:
    """Convert conformer to element coordinates

    :param symbols: Iterable of the symbols in the conformer
    :type symbols: Iterable[str]
    :param conformer: Conformer object
    :type conformer: Conformer
    :return: List of tuples with values element name, x, y and z coordinates
    :rtype: List[Tuple[str, float, float, float]]
    """
    result = []

    for atom, symbol in enumerate(symbols):
        p = conformer.GetAtomPosition(atom)
        result.append((symbol, p.x, p.y, p.z))

    return result


def _mol_to_coordinates(mol: "Mol", num_confs: int=10) -> List[Tuple[str, float, float, float]]:
    """Convert molecule object to list of coordinates

    :param mol: RDKit molecule object
    :type mol: Mol
    :param num_confs: Number of molecular conformers to generate, defaults to 10
    :type num_confs: int, optional
    :return: List of tuples of element name and x, y, z coordinates
    :rtype: List[Tuple[str, float, float, float]]
    """
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    conformer = get_conformer(mol=mol, num_confs=num_confs)

    return _conformer_to_coordinates(symbols=symbols, conformer=conformer)


def mol_to_xyz(mol: "Mol", num_confs: int=10) -> str:
    """Convert molecule object to XYZ file formatted string.

    :param mol: RDKit molecule object
    :type mol: Mol
    :param num_confs: Number of molecular conformers to generate, defaults to 10
    :type num_confs: int, optional
    :return: XYZ file formatted string
    :rtype: str
    """
    coordinates = _mol_to_coordinates(mol=mol, num_confs=num_confs)
    number_of_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)
    return coordinates_to_xyz(
        number_of_atoms=number_of_atoms,
        charge=charge,
        coordinates=coordinates
    )
