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
    terms : Dict = None
    system_indices : List[int] = None

    def __init__(self, data):
        self.__dict__ = data

    def add_terms(self, terms : List[Tuple[List[int], float]]) -> None:
        """ 
        Adds terms to the fermion hamiltonian.
        """
        logger.info(f"Adding {len(terms)} terms to fermion hamiltonian.")
        args = { 'hamiltonian': self.__dict__, 'fermion_terms': terms }
        args_json = json.dumps(map_tuples(args))
        result = qsharp.client._execute(f'%chemistry.fh.add_terms {args_json}', raise_on_stderr=True)
        self.__dict__ = result

class ProblemDescription(object):
    """
    Represents an electronic structure problem.
    """
    metadata : Dict  = None
    basis_set : Dict = None
    geometry : Dict = None
    coulomb_repulsion : Dict = None
    scf_energy : Dict = None
    scf_energy_offset : Dict = None
    fci_energy : Dict = None
    n_orbitals : int = None
    n_electrons : int = None
    energy_offset : Dict = None
    hamiltonian : Dict = None
    initial_state_suggestions : Dict = None

    def __init__ (self, data : Dict):
        self.__dict__ = data
    
    def load_fermion_hamiltonian(self, index_convention = 'UpDown') -> FermionHamiltonian:
        args = { 'problem_description': self.__dict__, 'index_convention': index_convention }
        args_json = json.dumps(map_tuples(args))
        data = qsharp.client._execute(f'%chemistry.fh.load {args_json}', raise_on_stderr=True)
        return FermionHamiltonian(data)

class Broombridge(object):
    """
    Represents an instance of a broombridge schema data
    """
    format : Dict = None
    generator : Dict = None
    bibliography : List[Dict] = None
    problem_description : List[ProblemDescription] = None

    def __init__(self, data: Dict):
        self.__dict__ = data
        # translate problem_descriptions into the actual datastructure:
        self.problem_description = [ ProblemDescription(p) for p in data["problem_description"] ]


def load_fermion_hamiltonian(filename: str) -> FermionHamiltonian:
    args = { 'file_name': filename }
    args_json = json.dumps(map_tuples(args))
    data = qsharp.client._execute(f'%chemistry.fh.load {args_json}', raise_on_stderr=True)
    return FermionHamiltonian(data)


def load_broombridge(filename: str) -> Broombridge:
    data = qsharp.client._execute(f'%chemistry.broombridge {filename}', raise_on_stderr=True)
    return Broombridge(data)

