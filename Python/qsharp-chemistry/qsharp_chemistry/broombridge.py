# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
"""This module specifies a class for representing a Boombridge data structure.
"""

import json
from typing import Dict

from .problem_description import ProblemDescription
from .fermion_hamiltonian import FermionHamiltonian

def _truncate(value: str):
    if len(value) > 100:
        return f"{value[:100]}..."
    return value

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

    def __repr__(self):
        return json.dumps(
            {
                k: v if k != "problem_description" else [
                    _truncate(str(_v)) for _v in v
                ] for k, v in self.__dict__.items()
            }, indent=4
        )
