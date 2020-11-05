# ---------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# ---------------------------------------------------------
"""This module specifies a class for representing a Boombridge data structure.
"""

from typing import Dict
from .problem_description import ProblemDescription
from .fermion_hamiltonian import FermionHamiltonian

class Broombridge(object):
    """
    Represents an instance of a Broombridge data structure
    """
    def __init__(self, data: Dict):
        self.__dict__ = data
        # translate problem_descriptions into the actual datastructure:
        self.problem_description = [ProblemDescription(p) for p in data["problem_description"]]

    def __eq__(self, other):
        if not isinstance(other, Broombridge):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.__dict__ == other.__dict__
