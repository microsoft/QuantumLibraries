# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for generating NWChem input deck format.
"""

import warnings
from typing import Union, TYPE_CHECKING

from qdk_chemistry.geometry import Geometry, format_geometry, format_geometry_from_mol
from qdk_chemistry.solvers.util import formatted_geometry_str, formatted_num_active_el

if TYPE_CHECKING:
    from rdkit.Chem.AllChem import Mol

__all__ = [
    "geometry_from_xyz",
    "create_input_deck"
]

_GEOMETRY_UNITS = ["au", "angstrom"]

# Template for generating an NWChem input deck
NW_CHEM_TEMPLATE = """
start {name}

echo
{memory}

geometry units {geometry_units}
symmetry c1
{geometry}
end

basis
* library {basis}
end

charge {charge}

scf
thresh {scf_thresh:.1e}
tol2e {scf_tol2e}
{rhf}
{spin}
{nopen}
end

tce
{method}
tilesize 1
2eorb
2emet 13
nroots {num_tce_root}
thresh {tce_thresh:.1e}
end

set tce:print_integrals T
set tce:qorb {num_orb}
set tce:qela {num_el_a}
set tce:qelb {num_el_b}

task tce {driver}
"""


def create_input_deck(
    mol: "Mol",
    mol_name: str, 
    num_active_orbitals: int,
    geometry: Union[str, Geometry] = None,
    memory: str = "memory stack 1000 mb heap 100 mb global 1000 mb noverify",
    geometry_units: str = "au",
    basis: str = "sto-3g",
    charge: int = 0,
    scf_thresh: float = 1.0e-10,
    scf_tol2e: float = 1.0e-10,
    rhf: str = "rhf",
    spin: str = "singlet",
    nopen: int = None,
    method: str = "ccsd",
    num_tce_root: int = 5,
    tce_thresh: float = 1.0e-6,
    driver: str = "energy",
    num_active_el: int = None,
) -> str:
    """Generate an NWChem input deck

    :param mol: RDKit Molecule object to use for calculating number of electrons if unspecified, defaults to None
    :type mol: Mol, optional
    :param mol_name: Molecule name
    :type mol_name: str
    :param geometry: Molecule geometry in the following format (each atom on a new line) : [Atom] [X] [Y] [Z]
    :type geometry: str
    :param num_active_orbitals: Number of orbitals in molecular ground state to set active space
    :type num_active_orbitals: int
    :param memory: Memory specification, defaults to "memory stack 1000 mb heap 100 mb global 1000 mb noverify"
    :type memory: str, optional
    :param geometry_units: Units used for geometry ("au"|"angstrom"), defaults to "au"
    :type geometry_units: str, optional
    :param basis: Basis to use for atoms, defaults to "sto-3g"
    :type basis: str, optional
    :param charge: Molecule charge, defaults to 0
    :type charge: int, optional
    :param scf_thresh: Threshold for SCF solver with one-decimal precision, defaults to 1.0e-10
    :type scf_thresh: float, optional
    :param scf_tol2e: 2-electron tolerance for SCF solver, defaults to 1.0e-10
    :type scf_tol2e: float, optional
    :param rhf: Restricted (Open-shell) Hartree Fock method. Either "rhf" or "rohf", defaults to "rhf"
    :type rhf: str, optional
    :param spin: Spin property, defaults to "singlet"
    :type spin: str, optional
    :param nopen: Number of singly occupied electron orbitals, defaults to None
    :type nopen: int, optional
    :param method: Calculation method, either "HF", "DFT" or "ccsd", defaults to "ccsd"
    :type method: str, optional
    :param num_tce_root: Number of excited states, defaults to 5
    :type num_tce_root: int, optional
    :param tce_thresh: Threshold for TCE solver with one-decimal precision, defaults to 1.0e-6
    :type tce_thresh: float, optional
    :param driver: Driver method, defaults to "energy"
    :type driver: str, optional
    :param num_active_el: Number of active electrons in molecule. This value is calculated based on atomic numbers if mol is provided. Defaults to None, defaults to None
    :type num_active_el: int, optional
    :raises ValueError: If neither Mol or num_active_el|geometry are specified.
    :return: NWChem input deck formatted string
    :rtype: str
    """
    assert geometry_units in _GEOMETRY_UNITS, f"Unknown geometry unit: {geometry_units}"

    num_active_el = formatted_num_active_el(mol=mol, num_active_el=num_active_el)
    geometry = formatted_geometry_str(mol=mol, geometry=geometry)
    nopen_str = f"nopen {nopen}" if nopen is not None else ""

    nw_chem = NW_CHEM_TEMPLATE.format(
        name=f"{mol_name}",
        memory=memory,
        geometry_units=geometry_units,
        geometry=geometry,
        basis=basis,
        charge=charge,
        scf_thresh=scf_thresh,
        scf_tol2e=scf_tol2e,
        rhf=rhf,
        spin=spin,
        nopen=nopen_str,
        method=method,
        num_tce_root=num_tce_root,
        tce_thresh=tce_thresh,
        num_orb=num_active_orbitals,
        num_el_a=num_active_el//2,
        num_el_b=num_active_el//2,
        driver=driver
    )

    return nw_chem
