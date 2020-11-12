from typing import Union

from qdk_chemistry.geometry import Geometry


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

{method}
"""

PSI4_TEMPLATE_MISC = """
{driver}('{method}')
"""

PSI4_TEMPLATE_FCIDUMP="""
e, wfn = {driver}('{method}', return_wfn=True)
fcidump(wfn, fname='fcidump', oe_ints=['EIGENVALUES'])
clean()
"""


def create_input_deck(
    mol: "Mol",
    mol_name: str = "",
    geometry: Union[str, Geometry] = None,
    charge: int = 0,
    spin: str = 1,
    basis: str = "ANO-RCC-MB",
    symmetry: str = "C1",
    integral_keyword: str = "",
    method: str = "SCF",
    driver: str = "energy",
    scf_type: str = "pk",
    memory_in_gb: int = 1,
    reference: str = "rhf",
    num_active_el: int = None,
    num_active_orbitals: int = 0,
    ci_root: int = 1,
    nsweeps: int = 5,
    max_bond_dimension: int = 500,
    d_convergence: float = 1e-8,
    e_convergence: float = 1e-8
) -> str:
    template = PSI4_TEMPLATE_FCIDUMP if driver == "energy" else PSI4_TEMPLATE_MISC
    method = template.format(driver=driver, method=method)

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
        method=method
    )
