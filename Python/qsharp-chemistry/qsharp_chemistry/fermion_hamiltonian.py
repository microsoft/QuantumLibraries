# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
"""This module contains the data structure for a fermion Hamiltonian as input to 
quantum chemistry algortihms.
"""

from dataclasses import dataclass
from typing import List, Tuple, Dict, Iterable

import logging
logger = logging.getLogger(__name__)

HTerm = Tuple[str, List[tuple]]

@dataclass
class FermionHamiltonian(object):
    """
    Represents a fermion Hamiltonian.
    """
    system_indices: List[int]
    terms: List[HTerm]

    @classmethod
    def from_dict(cls, data: dict):
        """Create FermionHamiltonian object from dictionary

        :param data: FermionHamiltonian dictionary
        :type data: dict
        :return: FermionHamiltonian object
        :rtype: FermionHamiltonian
        """
        return cls(**data)

    def __eq__(self, other):
        if not isinstance(other, FermionHamiltonian):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__

    def add_terms(self, fermion_terms: Iterable[Tuple[List[int], float]]) -> None:
        """
        Adds terms to the fermion Hamiltonian.
        """
        import qsharp
        logger.info(f"Adding {len(fermion_terms)} terms to fermion Hamiltonian.")
        result = qsharp.client._execute_magic(
            'chemistry.fh.add_terms',
            raise_on_stderr=True,
            hamiltonian=self.__dict__,
            fermion_terms=fermion_terms
        )
        self.__dict__ = result
