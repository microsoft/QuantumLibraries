from typing import Dict
from .problem_description import ProblemDescription
from .fermion_hamiltonian import FermionHamiltonian

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
