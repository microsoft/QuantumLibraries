# ---------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# ---------------------------------------------------------
"""This module contains the data structures for a quantum chemistry problem description.
This includes the fermion Hamiltonian and input state.
"""

import qsharp
import json
import typing

from typing import List, Tuple, Dict, Iterable
from enum import Enum

from .fermion_hamiltonian import FermionHamiltonian

import logging
logger = logging.getLogger(__name__)

class IndexConvention(Enum):
    UpDown = 1
    HalfUp = 2


class InputState(object):
    """
    Represents an input state.
    """
    def __init__(self, data: Dict):
        self.__dict__ = data

    def __eq__(self, other):
        if not isinstance(other, InputState):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__


class ProblemDescription(object):
    """
    Represents an electronic structure problem.
    """
    def __init__ (self, data: Dict):
        self.__dict__ = data

    def __eq__(self, other):
        if not isinstance(other, ProblemDescription):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__

    def load_fermion_hamiltonian(self, index_convention: IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
        """
        Loads the fermion Hamiltonian associated with this electronic structure problem.

        `index_convention` can be `UpDown` or `HalfUp`
        """
        logger.info(f"Loading fermion Hamiltonian from problem description using index_convention '{index_convention.name}'.")
        data = qsharp.client._execute_magic(
            'chemistry.fh.load',
            raise_on_stderr=True,
            problem_description=self.__dict__,
            index_convention=index_convention.name
        )
        return FermionHamiltonian(data)

    def load_input_state(self, wavefunction_label: str = '', index_convention: IndexConvention = IndexConvention.UpDown) -> FermionHamiltonian:
        """
        Loads the input state associated from this electronic structure problem with the corresponding label.

        If `wavefunction_label` is not specified, it loads the greedy (Hartree-Fock) state.

        `index_convention` can be 'UpDown' or 'HalfUp'
        """
        logger.info(f"Loading input state '{wavefunction_label}' from problem description using index_convention '{index_convention.name}'.")
        data = qsharp.client._execute_magic(
            'chemistry.inputstate.load',
            raise_on_stderr=True,
            problem_description=self.__dict__,
            wavefunction_label=wavefunction_label,
            index_convention=index_convention.name
        )
        return InputState(data)
