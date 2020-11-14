# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
"""This module contains shims for C#-defined functions for loading and encoding data via 
Broombridge-formatted files to solve quantum chemistry problems using quantum algorithms 
and are called via IQ# magics.
"""

import json
import typing

from typing import List, Tuple, Dict, Iterable
from enum import Enum

import qsharp

from .broombridge import Broombridge
from .fermion_hamiltonian import FermionHamiltonian
from .problem_description import IndexConvention, InputState

import logging
logger = logging.getLogger(__name__)

_CHEMISTRY_PACKAGE_NAME = "Microsoft.Quantum.Chemistry.Jupyter"

NumQubits = int
HamiltonianTermList = Tuple[List[Tuple[List[int], List[float]]]]
InputStateTerms = Tuple[int, List[Tuple[Tuple[float, float], List[int]]]]
EnergyOffset = float
JWEncodedData = Tuple[
    NumQubits,
    HamiltonianTermList,
    InputStateTerms,
    EnergyOffset
]

def enable_magic():
    """
    Enable the %chemistry magic command.
    """
    if _CHEMISTRY_PACKAGE_NAME not in [name for name, _ in list(qsharp.packages)]:
        qsharp.packages.add(_CHEMISTRY_PACKAGE_NAME)


def load_fermion_hamiltonian(file_name: str, index_convention: IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
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
    return FermionHamiltonian.from_dict(data)


def load_input_state(file_name: str, wavefunction_label: str = None, index_convention: IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
    """
    Loads the input state associated with the given label from the given file that contains Broombridge data.

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
    logger.info(f"Loading Broombridge data from '{file_name}'.")
    # NB: we don't use the execute_magic method here, since we don't need to
    #     JSON serialize any arguments in this case.
    data = qsharp.client._execute(f'%chemistry.broombridge {file_name}', raise_on_stderr=True)
    return Broombridge.from_dict(data)


def encode(hamiltonian: FermionHamiltonian, input_state: InputState) -> Tuple:
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


def load_and_encode(
    file_name: str,
    problem_description_index: int = 0,
    initial_state_label: str = None
) -> JWEncodedData:
    """Wrapper function for loading and encoding Broombridge file into
    JWEncodedData-compatible format.

    :param file_name: Broombridge file name
    :type file_name: str
    :param problem_description_index: Index of problem description to use, defaults to 0
    :type problem_description_index: int, optional
    :param initial_state_label: Label of initial state to use, defaults to first available label
    :type initial_state_label: str, optional
    """
    # TODO: This function accesses a file multiple times, find a way to cache 
    # or pass data explicitly
    # TODO: Make these tuples of lists etc. a bit more insightful data types
    broombridge_data =  load_broombridge(file_name)
    problem = broombridge_data.problem_description[problem_description_index]

    if initial_state_label is None:
        # Pick first in list
        initial_state_label = problem.initial_state_suggestions[0].get("Label")
        logger.info(f"Using initial state label: {initial_state_label}")

    input_state = load_input_state(file_name, initial_state_label)
    ferm_hamiltonian = problem.load_fermion_hamiltonian()
    num_qubits, hamiltonian_term_list, input_state_terms, energy_offset = encode(ferm_hamiltonian, input_state)

    return (num_qubits, hamiltonian_term_list, input_state_terms, energy_offset)
