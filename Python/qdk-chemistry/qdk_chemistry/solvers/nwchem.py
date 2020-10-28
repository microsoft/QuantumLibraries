from ..convert import num_electrons

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit.Chem.AllChem import Mol

# Template for generating an NW Chem input deck
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
thresh 1.0e-10
tol2e 1.0e-10
{rhf}
{spin}
end

tce
{method}
tilesize 1
2eorb
2emet 13
nroots 5
thresh 1.0e-6
end

set tce:print_integrals T
set tce:qorb {num_orb}
set tce:qela {num_el_a}
set tce:qelb {num_el_b}

task tce energy
"""

def geometry_from_xyz(xyz: str):
    """Generate geometry portion of NW Chem file from XYZ data

    Args:
        xyz (str): XYZ file format
    """
    lines = xyz.split("\n")
    return "\n".join(lines[2:])

def create_input_deck(
        mol_name: str, 
        geometry: str, 
        num_orbitals: int,
        memory: str = "memory stack 1000 mb heap 100 mb global 1000 mb noverify",
        geometry_units: str = "au", 
        basis: str = "sto-3g",
        charge: int = 0,
        rhf: str = "rhf",
        spin: str = "singlet",
        method: str = "ccsd",
        driver: str = "energy",
        mol: "Mol" = None,
        num_el: int = None,
    ):
    """Generate an NW Chem input deck

    Args:
        mol_name (str): Molecule name
        geometry (str): Molecule geometry in the following format (each atom on a new line):
            [Atom] [X] [Y] [Z]
        num_orbitals (int): Number of orbitals in molecular ground state to set active space
        memory (str, optional): Memory specification. Defaults to "memory stack 1000 mb heap 100 mb global 1000 mb noverify".
        geometry_units (str, optional): Units used for geometry. Defaults to "au".
        basis (str, optional): Basis to use for atoms. Defaults to "sto-3g".
        charge (int, optional): Molecule charge. Defaults to 0.
        rhf (str, optional): Restricted Hartree Fock method. Either "rhf" or "rohf". Defaults to "rhf".
        spin (str, optional): Spin property. Defaults to "singlet".
        method (str, optional): Calculation method, either "HF", "DFT" or "ccsd". Defaults to "ccsd".
        driver (str, optional): Driver method. Defaults to "energy".
        mol (Mol, optional): RDKit Molecule object to use for calculating number of electrons if unspecified. Defaults to None.
        num_el (int, optional): Manually enter the number of electrons. Defaults to None.

    Returns:
        str: NW Chem input deck formatted string
    """
    if mol is not None and num_el is None:
        num_el = num_electrons(mol)
    elif mol is None and num_el is None:
        raise ValueError("Cannot proceed: please provide either a Mol object or specify number of electrons")

    nw_chem = NW_CHEM_TEMPLATE.format(
        name=f"{mol_name}_test",
        memory=memory,
        geometry_units=geometry_units,
        geometry=geometry,
        basis=basis,
        charge=charge,
        rhf=rhf,
        spin=spin,
        method=method,
        num_orb=num_orbitals,
        num_el_a=num_el//2,
        num_el_b=num_el//2,
        driver=driver
    )

    return nw_chem
