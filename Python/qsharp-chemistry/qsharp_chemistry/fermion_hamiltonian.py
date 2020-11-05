from typing import List, Tuple, Dict, Iterable

import logging
logger = logging.getLogger(__name__)

class FermionHamiltonian(object):
    """
    Represents a fermion Hamiltonian.
    """
    def __init__(self, data):
        self.__dict__ = data

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
