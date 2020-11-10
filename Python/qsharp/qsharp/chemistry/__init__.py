# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
"""This module imports the qsharp_chemistry package such that it can be imported via qsharp.chemistry.
"""

try:
    import qsharp_chemistry
except ImportError:
    raise ImportError("Missing dependency: please run `pip install qsharp-chemistry` to use the qsharp.chemistry subpackage.")
else:
    from qsharp_chemistry.shims import load_input_state, load_fermion_hamiltonian, load_broombridge, encode
    from qsharp_chemistry.broombridge import Broombridge
    from qsharp_chemistry.fermion_hamiltonian import FermionHamiltonian
    from qsharp_chemistry.problem_description import ProblemDescription, InputState, IndexConvention
