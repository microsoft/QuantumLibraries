import os
import numpy as np
import pytest
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdmolops

from qdk_chemistry.xyz2mol.util import (
    get_atoms, 
    get_mol, 
    generate_structure_from_smiles, 
    get_proto_mol,
    read_xyz_file
)

from qdk_chemistry.xyz2mol.xyz2mol import xyz2mol
from qdk_chemistry.xyz2mol.ac import AC2mol

__TEST_SMILES__ = [
    'C[C-](c1ccccc1)C',
    'C[C-](C)c1ccccc1',
    'C=C([O-])CC',
    'C=C([NH3+])CC',
    'CC(=O)[O-]',
    'C[N+](=O)[O-]',
    'CS(CC)(=O)=O',
    'CS([O-])(=O)=O',
    'C=C(C)CC',
    'CC(C)CC',
    'C=C(N)CC',
    'C=C(C)C=C',
    'C#CC=C',
    'c1ccccc1',
    'c1ccccc1c1ccccc1',
    '[NH3+]CS([O-])(=O)=O',
    'CC(NC)=O',
    '[O-]c1ccccc1',
    'O=C(C=C1)C=CC1=CCC([O-])=O',
    'C#CC#C',
    'Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1',
    '[C+](C)(C)CC[C-](C)(C)',
    'O=C(C=C1)C=CC1=CCC([O-])=O',
    '[O-]c1ccccc1',
    'Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1',
    'CC1C=CC2C(C=CC2(C)C)=CC=1',
    'CC1=CC=C(C=CC2)C2C=C1',
    'CC1=CC=C(C2=CC=CC=C2)C=C1',
    'C1(CC2=CC=CC=C2)=CC=CC=C1',
    '[O-]c1ccccc1[O-]',
    'C[N+](=O)[O-]',
    'N#CC(C#N)=CC=C1C=CC=CC(=C1)c1ccc(cc1)[N+](=O)[O-]',
    'CNC([O-])=C([NH+]=C/CC(O)=O)C'
]

__TEST_FILES__ = [
    ("ethane.xyz", 0, "CC"),
    # ("examples/acetate.xyz", -1, "CC(=O)[O-]"),
    # ("examples/chiral_stereo_test.xyz", 0, "C/C=C/[C@@H](C)F"),
    # ("examples/propylbenzene.xyz", -1, "C[C-](C)c1ccccc1"),
]


@pytest.mark.parametrize("smiles", __TEST_SMILES__)
def test_smiles_from_adjacent_matrix(smiles):

    charged_fragments = True
    quick = True

    # Cut apart the smiles
    mol = get_mol(smiles)
    atoms = get_atoms(mol)
    charge = Chem.GetFormalCharge(mol)
    adjacent_matrix = Chem.GetAdjacencyMatrix(mol)

    #
    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    # Define new molecule template from atoms
    new_mol = get_proto_mol(atoms)

        # reconstruct the molecule from adjacent matrix, atoms and total charge
    new_mols = AC2mol(new_mol, adjacent_matrix, atoms, charge, charged_fragments, quick)
    
    new_mol_smiles_list = []
    for new_mol in new_mols:
        new_mol = Chem.RemoveHs(new_mol)
        new_mol_smiles = Chem.MolToSmiles(new_mol)

        new_mol_smiles_list.append(new_mol_smiles)

    assert canonical_smiles in new_mol_smiles_list

    return

@pytest.mark.parametrize("smiles", __TEST_SMILES__)
def test_smiles_from_coord_vdw(smiles):

    # The answer
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

    # generate forcefield coordinates
    atoms, coordinates = generate_structure_from_smiles(smiles)

    # Generate molobj from atoms, charge and coordinates
    mols = xyz2mol(atoms, coordinates, charge=charge)

    smiles_list = []
    for mol in mols:
    # For this test, remove chira. clean and canonical
        Chem.Kekulize(mol)
        mol = Chem.RemoveHs(mol)
        Chem.RemoveStereochemistry(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

        # Please look away. A small hack that removes the explicit hydrogens
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol)
        smiles_list.append(smiles)

    assert canonical_smiles in smiles_list

    return


@pytest.mark.parametrize("smiles", __TEST_SMILES__)
def test_smiles_from_coord_huckel(smiles):

    # The answer
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

    # generate forcefield coordinates
    atoms, coordinates = generate_structure_from_smiles(smiles)

    # Generate molobj from atoms, charge and coordinates
    mols = xyz2mol(atoms, coordinates, charge=charge, use_huckel=True)

    smiles_list = []
    for mol in mols:
        # For this test, remove chira. clean and canonical
        Chem.Kekulize(mol)
        mol = Chem.RemoveHs(mol)
        Chem.RemoveStereochemistry(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

        # Please look away. A small hack that removes the explicit hydrogens
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol)
        smiles_list.append(smiles)

    assert canonical_smiles in smiles_list

    return


@pytest.mark.parametrize("filename, charge, answer", __TEST_FILES__)
def test_smiles_from_xyz_files(filename, charge, answer):

    charged_fragments = True
    quick = True

    atoms, charge_read, coordinates = read_xyz_file(os.path.join(os.path.split(__file__)[0], filename))
    mols = xyz2mol(atoms, coordinates, charge=charge)

    smiles_list = []
    for mol in mols:
        mol = Chem.RemoveHs(mol)

        smiles = Chem.MolToSmiles(mol)
        smiles_list.append(smiles)

    print(answer, smiles_list)

    assert answer in smiles_list

    return
