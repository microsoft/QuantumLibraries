from rdkit.Chem import AllChem as Chem
from qdk_chemistry.geometry import Element, Geometry

def test_Element():
    c1 = Element(
        name="C",
        x=0.0,
        y=0.0,
        z=0.0
    )

    c2 = Element(
        name="C",
        x=1.0,
        y=0.0,
        z=0.0
    )

    g = Geometry([c1, c2])
    assert g[0] == c1
    assert g[1] == c2


def test_geometry_from_mol():
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    g = Geometry.from_mol(mol)
    assert g.charge == 0
    assert len(g) == 3
    assert [el.name for el in g] == ["O", "H", "H"]

def test_geometry_from_xyz():
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    g = Geometry.from_mol(mol)
    xyz = g.to_xyz()
    gprime = Geometry.from_xyz(xyz)

    assert g == gprime
