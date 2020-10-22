import json
import typing

from typing import List, Tuple, Dict, Iterable
from enum import Enum

import qsharp

from .broombridge import Broombridge
from .fermion_hamiltonian import FermionHamiltonian, InputState
from .problem_description import IndexConvention

import logging
logger = logging.getLogger(__name__)

_CHEMISTRY_NAMESPACE = "Microsoft.Quantum.Chemistry.Jupyter"

def enable_magic():
    """
    Enable the %chemistry magic command.
    """
    if _CHEMISTRY_NAMESPACE not in [name for name, _ in list(qsharp.packages)]:
        qsharp.packages.add(_CHEMISTRY_NAMESPACE)


def load_fermion_hamiltonian(file_name: str, index_convention : IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
    """
    Loads the fermion Hamiltonian from the given file that contains broombridge data.

    `index_convention` can be 'UpDown' or 'HalfUp'
    """
    logger.info(f"Loading fermion Hamiltonian from '{file_name}' using index_convention '{index_convention.name}'.")
    enable_magic()
    data = qsharp.client._execute_magic(
        'chemistry.fh.load',
        raise_on_stderr=True,
        file_name=file_name,
        index_convention=index_convention.name
    )
    return FermionHamiltonian(data)


def load_input_state(file_name: str, wavefunction_label : str = None, index_convention : IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
    """
    Loads the input state associated with the given labe from the given file that contains broombridge data..

    If `wavefunction_label` is not specified, it loads the greedy (Hartree–Fock) state.

    `index_convention` can be 'UpDown' or 'HalfUp'
    """
    logger.info(f"Loading input state '{wavefunction_label}' from '{file_name}' using index_convention '{index_convention.name}'.")
    enable_magic()
    data = qsharp.client._execute_magic(
        'chemistry.inputstate.load',
        raise_on_stderr=True,
        file_name=file_name,
        wavefunction_label=wavefunction_label,
        index_convention=index_convention.name
    )
    return InputState(data)


def load_broombridge(file_name: str) -> Broombridge:
    """
    Loads a Broombridge data file.
    """
    enable_magic()
    logger.info(f"Loading broombridge data from '{file_name}'.")
    # NB: we don't use the execute_magic method here, since we don't need to
    #     JSON serialize any arguments in this case.
    data = qsharp.client._execute(f'%chemistry.broombridge {file_name}', raise_on_stderr=True)
    return Broombridge(data)


def encode(hamiltonian : FermionHamiltonian, input_state : InputState) -> Tuple:
    """
    Encodes the given Hamiltonian and input state using the Jordan Wigner encoding
    that can be used to run chemistry simulations using Q#'s chemistry library.
    """
    enable_magic()
    logger.info(f"Doing jw encoding.")
    data = qsharp.client._execute_magic(
        'chemistry.encode',
        raise_on_stderr=True,
        hamiltonian=hamiltonian.__dict__,
        input_state=input_state.__dict__
    )
    return data
