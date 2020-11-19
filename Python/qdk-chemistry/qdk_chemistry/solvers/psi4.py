import enum
from typing import Union

from qdk_chemistry.geometry import Geometry

class Basis(enum.Enum):
    gaussian = "3-21G"
    pople = "6-31G"
    dunning = "cc-pVDZ"
    karlsruhe = "def2-SVP"


PSI4_TEMPLATE = """
memory {N} GB

molecule {name} {{
symmetry {symmetry}
{charge} {spin}
  {geometry}
}}

set {{
  basis {basis}
  scf_type {scf_type}
  reference {reference}
  d_convergence 1e-8
  e_convergence 1e-8
}}

{method_section}
"""

PSI4_TEMPLATE_ENERGY = "e = {driver}('{method}')"


PSI4_TEMPLATE_FCIDUMP = """
e, wfn = {driver}('SCF', return_wfn=True)
{optional_energy}
fcidump(wfn, fname='fcidump', oe_ints=['EIGENVALUES'])
clean()
"""


def create_input_deck(
    mol: "Mol",
    mol_name: str = "",
    geometry: Union[str, Geometry] = None,
    charge: int = 0,
    spin: str = 1,
    basis: str = "3-21G",
    symmetry: str = "C1",
    method: str = "SCF",
    driver: str = "energy",
    scf_type: str = "PK",
    memory_in_gb: int = 1,
    reference: str = "rhf"
) -> str:
    """Create Psi4 input formatted string

    :param mol: RDKit Molecule object
    :type mol: Mol
    :param mol_name: Molecule name, defaults to ""
    :type mol_name: str, optional
    :param geometry: Geometry object, defaults to None
    :type geometry: Union[str, Geometry], optional
    :param charge: Molecule charge, defaults to 0
    :type charge: int, optional
    :param spin: Molecule spin, defaults to 1
    :type spin: str, optional
    :param basis: Molecule basis, defaults to "3-21G"
    :type basis: str, optional
    :param symmetry: Molecule symmetry (example options: "C1", "C2v", "D6h", "Cs"),
         defaults to "C1". When requesting Broombridge, make sure to set this to "C1"
         because Broombridge does not support group symmetry
    :type symmetry: str, optional
    :param method: Molecule method (options: "SCF", "HF", "DFT", "CC2", "CCSD", "CASSCF", 
        "MP2"), defaults to "SCF"
    :type method: str, optional
    :param driver: Driver method ("energy" or "optimize"), defaults to "energy"
    :type driver: str, optional
    :param scf_type: SCF solver type ("DIRECT", "DF", "PK", "OUT_OF_CORE" or "PS"), defaults to "PK"
    :type scf_type: str, optional
    :param memory_in_gb: Memory used in GB
    :type memory_in_gb: int
    :param reference: Reference, defaults to "rhf"
    :type reference: str, optional
    :return: Psi4 input-formatted string
    :rtype: str
    """
    assert basis in Basis._value2member_map_, f"Unknown basis: {basis}. Please choose from the following: {Basis._value2member_map_.values()}"
    optional_energy = PSI4_TEMPLATE_ENERGY.format(driver=driver, method=method) if method.upper() != "SCF" else ""
    template = PSI4_TEMPLATE_FCIDUMP.format(driver=driver, optional_energy=optional_energy)
    method_section = template.format(driver=driver, method=method)

    return PSI4_TEMPLATE.format(
        N=memory_in_gb,
        name=mol_name,
        symmetry=symmetry,
        charge=charge,
        spin=spin,
        geometry=geometry,
        basis=basis,
        scf_type=scf_type,
        reference=reference,
        method_section=method_section
    )
