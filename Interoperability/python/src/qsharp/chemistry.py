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

    def __init__(self, terms, indices):
        self.terms = terms
        self.system_indices = indices

    def add_terms(self, terms : List[Tuple[List[int], float]]) -> None:
        """ 
        Adds terms to the fermion hamiltonian.
        """
        logger.info(f"Adding {len(terms)} terms to fermion hamiltonian.")
        args = { 'hamiltonian': { 'terms': self.terms, 'system_indices': self.system_indices }, 'fermion_terms': terms }
        args_json = json.dumps(map_tuples(args))
        result = qsharp.client._execute(f'%chemistry.fh.add_terms {args_json}', raise_on_stderr=True)
        self.terms = result["terms"]
        self.system_indices = result["system_indices"]

class ProblemDescription(object):
    """
    Represents an electronic structure problem.
    """

    def __init__ (self, data : Dict):
        self.metadata  = data["metadata"]
        self.basis_set = data["basis_set"]
        self.geometry = data["geometry"]
        self.coulomb_repulsion = data["coulomb_repulsion"]
        self.scf_energy = data["scf_energy"]
        self.scf_energy_offset = data["scf_energy_offset"]
        self.fci_energy = data["fci_energy"]
        self.n_orbitals = data["n_orbitals"]
        self.n_electrons = data["n_electrons"]
        self.energy_offset = data["energy_offset"]
        self.hamiltonian = data["hamiltonian"]
        self.initial_state_suggestions = data["initial_state_suggestions"]
    

class Broombridge(object):
    """
    Represents an instance of a broombridge schema data
    """
    format = None
    generator = None
    bibliography = None
    problem_descriptions = None

    def __init__(self, format: Dict, generator: Dict, bibliography: List[Dict], problem_descriptions: List[Dict]):
        self.format = format
        self.generator = generator
        self.bibliography = bibliography
        self.problem_descriptions = problem_descriptions


def load_fermion_hamiltonian(filename: str):
    args = { 'file_name': filename }
    fh = qsharp.client._execute(f'%chemistry.fh.load {json.dumps(args)}', raise_on_stderr=True)
    
    return FermionHamiltonian(fh["terms"], fh["system_indices"])


def load_broombridge(filename: str):
    bb = qsharp.client._execute(f'%chemistry.fh.load {filename}', raise_on_stderr=True)    
    return Broombridge(bb["Format"], bb["Generator"], bb["Bibliography"], bb["ProblemDescriptions"])

