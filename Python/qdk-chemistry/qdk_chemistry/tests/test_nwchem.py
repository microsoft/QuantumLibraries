# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
import pytest
from unittest import mock

from qdk_chemistry.solvers.nwchem import create_input_deck, format_geometry_from_mol

from rdkit.Chem import AllChem as Chem


@pytest.fixture()
def test_deck():
    return """
start HHO_test

echo
memory stack 1000 mb heap 100 mb global 1000 mb noverify

geometry units au
symmetry c1
O 0.002 0.398 0.0
H 0.762 -0.203 0.0
H -0.764 -0.195 0.0
end

basis
* library sto-3g
end

charge 0

scf
thresh 1.0e-10
tol2e 1e-10
rhf
singlet

end

tce
ccsd
tilesize 1
2eorb
2emet 13
nroots 5
thresh 1.0e-06
end

set tce:print_integrals T
set tce:qorb 7
set tce:qela 5
set tce:qelb 5

task tce energy
"""


def test_nwchem(geometry, h2o, test_deck):
    mol_name = "HHO_test"
    with mock.patch("qdk_chemistry.geometry.Geometry.from_mol") as _m:
        _m.return_value = geometry

        nw_chem_input = create_input_deck(
            mol_name=mol_name,
            num_active_orbitals=7,
            mol=h2o
        )

    assert nw_chem_input == test_deck


def test_nwchem_pass_geometry(geometry, h2o, test_deck):
    mol_name = "HHO_test"
    nw_chem_input = create_input_deck(
        mol_name=mol_name,
        geometry=geometry,
        num_active_orbitals=7,
        mol=h2o
    )

    assert nw_chem_input == test_deck
