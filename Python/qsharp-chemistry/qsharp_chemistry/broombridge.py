# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
"""This module specifies a class for representing a Boombridge data structure.
"""

import json
from typing import Dict, List
from dataclasses import dataclass

from .problem_description import ProblemDescription
from .fermion_hamiltonian import FermionHamiltonian

def _truncate(value: str):
    if len(value) > 100:
        return f"{value[:100]}..."
    return value

@dataclass
class Broombridge(object):
    """
    Represents an instance of a Broombridge data structure
    """
    schema: str
    format: dict
    generator: dict
    bibliography: list
    problem_description: List[ProblemDescription]

    @classmethod
    def from_dict(cls, data: Dict):
        """Create Boombridge object from dictionary

        :param data: Dict representation of Broombridge object
        :type data: Dict
        :return: Broombridge object
        :rtype: Broombridge
        """
        problem_description = [ProblemDescription.from_dict(p) for p in data["problem_description"]]
        kwargs = {k.lstrip("$"): v for k, v in data.items() if k != problem_description}
        kwargs["problem_description"] = problem_description
        obj = cls(**kwargs)
        # translate problem_descriptions into the actual datastructure:
        return obj

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
