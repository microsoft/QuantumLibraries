# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for converting to and from RDKit data structures
"""

from typing import Iterable, TYPE_CHECKING

from rdkit.Chem import AllChem as Chem

if TYPE_CHECKING:
    from rdkit.Chem import Mol, Conformer

def conformer_to_xyz(number_of_atoms: int, charge: int, symbols: Iterable[str], conformer: "Conformer") -> str:
    """Convert conformer to XYZ file formatted string.

    :param number_of_atoms: Number of atoms in the conformer
    :type number_of_atoms: int
    :param charge: Charge of the conformer
    :type charge: int
    :param symbols: Iterable of the symbols in the conformer
    :type symbols: Iterable[str]
    :param conformer: Conformer object
    :type conformer: Conformer
    :return: XYZ file format
    :rtype: str
    """
    result = [
        f"{number_of_atoms}",
        "title"
    ]

    for atom, symbol in enumerate(symbols):
        p = conformer.GetAtomPosition(atom)
        result.append(f"{symbol} {str(p.x)} {str(p.y)} {str(p.z)} ")

    if charge != 0:
        result.extend([
            "$set",
            f"chrg {charge}",
            "$end"
        ])

    return "\n".join(result)

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

def mol_to_xyz(mol: "Mol", num_confs: int=10) -> str:
    """Convert molecule object to XYZ file formatted string.

    :param mol: RDKit molecule object
    :type mol: Mol
    :param num_confs: Number of molecular configurations to generate, defaults to 10
    :type num_confs: int, optional
    :return: XYZ file formatted string
    :rtype: str
    """
    number_of_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    conformer = get_conformer(mol=mol, num_confs=num_confs)

    return conformer_to_xyz(number_of_atoms, charge, symbols, conformer)


def num_electrons(mol: "Mol") -> int:
    """Calculate the number of electrons in the molecule

    :param mol: RDKit molecule object
    :type mol: Mol
    :return: Number of electrons or sum of atomic number of atoms in molecule
    :rtype: int
    """
    return sum([atom.GetAtomicNum() for atom in mol.GetAtoms()])
