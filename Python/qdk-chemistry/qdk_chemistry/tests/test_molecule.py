import os
from qdk_chemistry.molecule import Molecule

def test_molecule(tmp_path):
    path, _ = os.path.split(__file__)
    molecule = Molecule.from_xyz(os.path.join(path, "h2o.xyz"))
    molecule.create_input(
        molecule_name="H2O",
        file_name="h2o.nw",
        solver="NWChem",
        base_path=tmp_path,
        num_active_orbitals=2
    )

    assert os.path.isfile(os.path.join(tmp_path, "h2o.nw"))
