#!/bin/env python
# -*- coding: utf-8 -*-
##
# chemistry.py: enables using Q# chemistry library from Python.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

## IMPORTS ##

import qsharp
import json
import typing

from qsharp.serialization import map_tuples
from typing import List, Tuple, Dict

## LOGGING ##

import logging
logger = logging.getLogger(__name__)

## Enable the magic:
qsharp.packages.add("microsoft.quantum.chemistry.jupyter")


## CLASSES ##

class FermionHamiltonian(object):
    """
    Represents a fermion hamiltonian.
    """
    def __init__(self, data):
        self.__dict__ = data

    def __eq__(self, other): 
        if not isinstance(other, FermionHamiltonian):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__

    def add_terms(self, fermion_terms : List[Tuple[List[int], float]]) -> None:
        """ 
        Adds terms to the fermion hamiltonian.
        """
        logger.info(f"Adding {len(fermion_terms)} terms to fermion hamiltonian.")
        args = { 'hamiltonian': self.__dict__, 'fermion_terms': fermion_terms }
        args_json = json.dumps(map_tuples(args))
        result = qsharp.client._execute(f'%chemistry.fh.add_terms {args_json}', raise_on_stderr=True)
        self.__dict__ = result

class InputState(object):
    """
    Represents an input state.
    """
    def __init__(self, data: Dict):
        self.__dict__ = data

    def __eq__(self, other): 
        if not isinstance(other, FermionHamiltonian):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__

class ProblemDescription(object):
    """
    Represents an electronic structure problem.
    """
    def __init__ (self, data : Dict):
        self.__dict__ = data
    
    def __eq__(self, other): 
        if not isinstance(other, FermionHamiltonian):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__

    def load_fermion_hamiltonian(self, index_convention = 'UpDown') -> FermionHamiltonian:
        """
        Loads the fermion hamiltonian associated with this electronic structure problem.

        `index_convention` can be 'UpDown' or 'HalfUp'
        """
        logger.info(f"Loading fermion hamiltonian from problem description using index_convention '{index_convention}'.")
        args = { 'problem_description': self.__dict__, 'index_convention': index_convention }
        args_json = json.dumps(map_tuples(args))
        data = qsharp.client._execute(f'%chemistry.fh.load {args_json}', raise_on_stderr=True)
        return FermionHamiltonian(data)

        
    def load_input_state(self, wavefunction_label ='', index_convention = 'UpDown') -> FermionHamiltonian:
        """
        Loads the input state associated from this electronic structure problem with the corresponding label.

        If `wavefunction_label` is not specified, it loads the greedy (Hartree-Fock) state.

        `index_convention` can be 'UpDown' or 'HalfUp'
        """
        logger.info(f"Loading input state '{wavefunction_label}' from problem description using index_convention '{index_convention}'.")
        args = { 'problem_description': self.__dict__, 'wavefunction_label': wavefunction_label, 'index_convention': index_convention }
        args_json = json.dumps(map_tuples(args))
        data = qsharp.client._execute(f'%chemistry.inputstate.load {args_json}', raise_on_stderr=True)
        return InputState(data)


class Broombridge(object):
    """
    Represents an instance of a broombridge schema data
    """
    def __init__(self, data: Dict):
        self.__dict__ = data
        # translate problem_descriptions into the actual datastructure:
        self.problem_description = [ ProblemDescription(p) for p in data["problem_description"] ]

    def __eq__(self, other): 
        if not isinstance(other, FermionHamiltonian):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__


def load_fermion_hamiltonian(file_name: str, index_convention = 'UpDown') -> FermionHamiltonian:
    logger.info(f"Loading fermion hamiltonian from '{file_name}' using index_convention '{index_convention}'.")
    """
    Loads the fermion hamiltonian from the given file that contains broombridge data.

    `index_convention` can be 'UpDown' or 'HalfUp'
    """
    args = { 'file_name': file_name, 'index_convention': index_convention }
    args_json = json.dumps(map_tuples(args))
    data = qsharp.client._execute(f'%chemistry.fh.load {args_json}', raise_on_stderr=True)
    return FermionHamiltonian(data)


def load_input_state(file_name: str, wavefunction_label : str, index_convention = 'UpDown') -> FermionHamiltonian:
    """
    Loads the input state associated with the given labe from the given file that contains broombridge data..

    If `wavefunction_label` is not specified, it loads the greedy (Hartree-Fock) state.

    `index_convention` can be 'UpDown' or 'HalfUp'
    """
    logger.info(f"Loading input state '{wavefunction_label}' from '{file_name}' using index_convention '{index_convention}'.")
    args = { 'file_name': file_name, 'wavefunction_label': wavefunction_label, 'index_convention': index_convention }
    args_json = json.dumps(map_tuples(args))
    data = qsharp.client._execute(f'%chemistry.inputstate.load {args_json}', raise_on_stderr=True)
    return FermionHamiltonian(data)


def load_broombridge(file_name: str) -> Broombridge:
    """
    Loads a broombridge data file.
    """
    logger.info(f"Loading broombridge data from '{file_name}'.")
    data = qsharp.client._execute(f'%chemistry.broombridge {file_name}', raise_on_stderr=True)
    return Broombridge(data)

    
def encode(hamiltonian : FermionHamiltonian, input_state : InputState) -> Tuple:
    """
    Encodes the given hamiltonian and input state using the Jordan Wigner encoding
    that can be used to run chemistry simulations using Q#'s chemistry library.
    """
    logger.info(f"Doing jw encoding.")
    args = { 'hamiltonian': hamiltonian.__dict__, 'input_state': input_state.__dict__ }
    args_json = json.dumps(map_tuples(args))
    data = qsharp.client._execute(f'%chemistry.encode {args_json}', raise_on_stderr=True)
    return data

