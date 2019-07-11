#!/bin/env python
# -*- coding: utf-8 -*-
##
# types.py: Provides types for interoperability with Q#.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

from enum import IntEnum
import random

try:
    import qutip as qt
except ImportError:
    qt = None

## ENUMS ######################################################################

class Result(IntEnum):
    """
    Represents the `Result Q# type <https://docs.microsoft.com/quantum/language/type-model#primitive-types>`_.
    """
    Zero = 0
    One = 1

    @classmethod
    def sample(cls):
        """
        Samples a value of type Result uniformly at random.
        """
        return random.choice(list(cls))

class Pauli(IntEnum):
    """
    Represents the `Pauli Q# type <https://docs.microsoft.com/quantum/language/type-model#primitive-types>`_.
    """
    I = 0
    X = 1
    Y = 2
    Z = 3

    @classmethod
    def sample(cls):
        """
        Samples a value of type Pauli uniformly at random.
        """
        return random.choice(list(cls))

    if qt is not None:
        def as_qobj(self) -> qt.Qobj:
            """
            Represents this value as a QuTiP Qobj (quantum object).

            :return: A representation of this Pauli value as a QuTiP quantum
                 object.
            """
            if self == Pauli.I:
                return qt.qeye(2)
            elif self == Pauli.X:
                return qt.sigmax()
            elif self == Pauli.Y:
                return qt.sigmay()
            elif self == Pauli.Z:
                return qt.sigmaz()
            else:
                raise ValueError(f"Unrecognized Pauli value {self}.")
