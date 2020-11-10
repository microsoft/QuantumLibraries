# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
import pytest

from rdkit.Chem import AllChem as Chem

from qdk_chemistry.geometry import Geometry, Element


@pytest.fixture()
def h2o():
    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    return mol


@pytest.fixture()
def geometry(h2o):
    el = [
        Element("O", 0.002, 0.398, 0.0),
        Element("H", 0.762, -0.203, 0.0),
        Element("H", -0.764, -0.195, 0.0)
    ]
    return Geometry(el, charge=Chem.GetFormalCharge(h2o))
