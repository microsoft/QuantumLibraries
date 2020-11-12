# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
import pytest
from unittest import mock

from qdk_chemistry.solvers.psi4 import create_input_deck
from qdk_chemistry.geometry import format_geometry_from_mol

from rdkit.Chem import AllChem as Chem


@pytest.fixture()
def test_deck():
    return """
memory 1 GB

molecule  {
symmetry C1
0 1
  None
}

set {
  basis ANO-RCC-MB
  scf_type pk
  reference rhf
  d_convergence 1.0e-08
  e_convergence 1.0e-08
}


e, wfn = energy('SCF', return_wfn=True)
fcidump(wfn, fname='fcidump', oe_ints=['EIGENVALUES'])
clean()

"""


def test_psi4(geometry, h2o, test_deck):
    mol_name = "HHO_test"
    with mock.patch("qdk_chemistry.geometry.Geometry.from_mol") as _m:
        _m.return_value = geometry

        psi4_input = create_input_deck(
            mol=h2o,
            method="SCF",
            driver="energy"
        )

    assert psi4_input == test_deck
