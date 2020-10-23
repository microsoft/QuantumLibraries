from typing import Iterable, TYPE_CHECKING

from rdkit.Chem import AllChem as Chem

if TYPE_CHECKING:
    from rdkit.Chem import Mol, Conformer

def conformer_to_xyz(number_of_atoms: int, charge: int, symbols: Iterable[str], conformer: "Conformer") -> str:
    """Convert conformer to XYZ file formatted string.

    Args:
        number_of_atoms (int): Number of atoms in the conformer
        charge (int): Charge of the conformer
        symbols (Iterable[str]): Iterable of the symbols in the conformer
        conformer (Conformer): Conformer object

    Returns:
        str: XYZ file format
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

    Args:
        mol (Mol): Molecular fragment
        num_confs (int, optional): [description]. Defaults to 10.

    Returns:
        "Conformer": Lowest-energy conformer for molecule
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

    Args:
        mol (Mol): [description]
        num_confs (int, optional): [description]. Defaults to 10.

    Returns:
        str: XYZ file formatted string
    """
    number_of_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    conformer = get_conformer(mol=mol, num_confs=num_confs)

    return conformer_to_xyz(number_of_atoms, charge, symbols, conformer)
