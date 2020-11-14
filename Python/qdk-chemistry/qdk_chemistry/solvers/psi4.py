from typing import Union

from qdk_chemistry.geometry import Geometry

class Basis():
    pople = "6-31G"
    dunning = "cc-pVDZ",
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
  d_convergence {d_convergence:.1e}
  e_convergence {e_convergence:.1e}
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
    driver: str = "energy", # "energy", "optimize"
    scf_type: str = "PK", # "DIRECT", "DF", "PK", "OUT_OF_CORE", "PS"
    memory_in_gb: int = 1,
    reference: str = "rhf",
    d_convergence: float = 1e-8,
    e_convergence: float = 1e-8
) -> str:
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
        d_convergence=d_convergence,
        e_convergence=e_convergence,
        method_section=method_section
    )
