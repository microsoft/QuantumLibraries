# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
import pytest
from unittest import mock

from qdk_chemistry.solvers.openmolcas import create_input_deck
from qdk_chemistry.geometry import format_geometry_from_mol

from rdkit.Chem import AllChem as Chem


@pytest.fixture()
def test_deck():
    return """
&GATEWAY
Coord
  3
  
  O 0.002 0.398 0.0
  H 0.762 -0.203 0.0
  H -0.764 -0.195 0.0
Basis=ANO-RCC-MB
Group=C1

&SEWARD


&SCF
Charge=0
Spin=1

&RASSCF
  DMRG
  FCIDUMP
  TYPEINDEX
  Spin=1
  Charge=0

  RGINPUT
  nsweeps = 5
  max_bond_dimension = 500
  ENDRG

"""


def test_openmolcas(geometry, h2o, test_deck):
    mol_name = "HHO_test"
    with mock.patch("qdk_chemistry.geometry.Geometry.from_mol") as _m:
        _m.return_value = geometry

        openmolcas_input = create_input_deck(
            charge=0,
            spin=1,
            method="fcidump",
            mol=h2o
        )

    assert openmolcas_input == test_deck
