# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for generating OpenMolcas input deck format.
"""

from typing import Union

from qdk_chemistry.geometry import Geometry
from qdk_chemistry.solvers.util import (
    formatted_geometry_str, 
    formatted_num_active_el, 
    num_atoms_from_mol
)

from rdkit.Chem.rdmolops import GetFormalCharge


GEOMETRY_LINE_SEP = "\n  "

OPENMOLCAS_TEMPLATE = """
&GATEWAY
Coord
  {num_atoms}
  {mol_name}
  {geometry}
Basis={basis}
Group={symmetry}

&SEWARD
{integral_keyword}

&SCF
Charge={charge}
Spin={spin}
{rasscf}
"""

OPENMOLCAS_TEMPLATE_CASSCF = """
&RASSCF
  Charge= {charge}
  Nactel  =  {num_active_el}
  Ras2  =  {num_active_orbitals}
  Ciroot  =  {nroot} {nroot} 1
"""

OPEN_MOLCAS_TEMPLATE_FCIDUMP = """
&RASSCF
  DMRG
  FCIDUMP
  TYPEINDEX
  Spin={spin}
  Charge={charge}

  RGINPUT
  nsweeps = {nsweeps}
  max_bond_dimension = {max_bond_dimension}
  ENDRG
"""


def create_input_deck(
    mol: "Mol",
    mol_name: str = "",
    geometry: Union[str, Geometry] = None,
    charge: int = None,
    spin: str = 1,
    basis: str = "ANO-RCC-MB",
    symmetry: str = "C1",
    integral_keyword: str = "",
    method: str = "",
    num_active_el: int = None,
    num_active_orbitals: int = 0,
    ci_root: int = 1,
    nsweeps: int = 5,
    max_bond_dimension: int = 500
) -> str:
    """Convenience function for creating an OpenMolcas input formatted string.

    :param mol: Molecule object
    :type mol: Mol
    :param mol_name: Molecule name, defaults to ""
    :type mol_name: str, optional
    :param geometry: Molecule geometry, defaults to None
    :type geometry: Union[str, Geometry], optional
    :param charge: Charge, defaults to None
    :type charge: int, optional
    :param spin: Spin, defaults to 1
    :type spin: str, optional
    :param basis: Molecule basis, defaults to "ANO-RCC-MB"
    :type basis: str, optional
    :param symmetry: Molecule summetry, defaults to "C1"
    :type symmetry: str, optional
    :param integral_keyword: Integral keyword, e.g. "Cholesky", defaults to ""
    :type integral_keyword: str, optional
    :param method: Method to use, e.g. "FCIDUMP" or "CASSCF", defaults to ""
    :type method: str, optional
    :param num_active_el: Number of active electrons, defaults to None (will be calculated from molecule if unspecified)
    :type num_active_el: int, optional
    :param num_active_orbitals: number of orbitals in each symmetry for the RAS2 orbital subspace, defaults to 0
    :type num_active_orbitals: int, optional
    :param ci_root: CI root(s), see https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rasscf.html#optional-keywords, defaults to 1 (ground state)
    :type ci_root: int, optional
    :param nsweeps: Number of sweeps for FCIDUMP method, defaults to 5
    :type nsweeps: int, optional
    :param max_bond_dimension: Maximum bond dimension for FCIDUMP method, defaults to 500
    :type max_bond_dimension: int, optional
    :return: OpenMolcas input string
    :rtype: str
    """
    formatted_geometry = formatted_geometry_str(
        mol=mol, 
        geometry=geometry, 
        line_sep=GEOMETRY_LINE_SEP
    )
    num_atoms = len(geometry) if isinstance(geometry, Geometry) else num_atoms_from_mol(mol)
    charge = charge if charge is not None else GetFormalCharge(mol)

    if method.upper() == "CASSCF":
        num_active_el = formatted_num_active_el(mol=mol, num_active_el=num_active_el)
        rasscf = OPENMOLCAS_TEMPLATE_CASSCF.format(
            charge=charge,
            num_active_el=num_active_el,
            num_active_orbitals=num_active_orbitals,
            nroot=ci_root
        )

    elif method.upper() == "FCIDUMP":
        rasscf = OPEN_MOLCAS_TEMPLATE_FCIDUMP.format(
            spin=spin,
            charge=charge,
            nsweeps=nsweeps,
            max_bond_dimension=max_bond_dimension
        )

    else:
      rasscf = ""

    return OPENMOLCAS_TEMPLATE.format(
        num_atoms=num_atoms,
        mol_name=mol_name,
        geometry=formatted_geometry,
        basis=basis,
        symmetry=symmetry,
        integral_keyword=integral_keyword,
        charge=charge,
        spin=spin,
        rasscf=rasscf
    )
