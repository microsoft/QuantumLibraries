[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

# QAOA in Q#

The project is still in progress. This readme will be extended as the project develops.

This project provides a hybrid quantum-classical algorithm for solving optimization problems. 
It includes a Q# implementation of the Quantum Approximate Optimization Algorithm ([QAOA](https://arxiv.org/abs/1411.4028)) together with a classical optimizer in C#.
The classical optimizer uses a quantum objective function to choose hopefully optimal parameters for running the QAOA.
The quantum objective function is currently evaluated on a simulator backend provided by Q#.

How to run it?
1) Import this project to Microsoft Visual Studio or similar.
2) Use Driver.cs to prepare your ProblemInstance and run the project (some examples also provided).

Current limitations:

- an optimization problem shall be encoded into a Hamiltonian consisting of Z operators,
- support for up to 2-local Hamiltonians,
- input consists of arrays of coefficients for 1-local and 2-local Hamiltonian terms,
- a gradient-free Cobyla optimizer is used for finding good QAOA input parameters.

The high-level diagram of the implementation (notation comes from the [QAOA paper](https://arxiv.org/abs/1411.4028)):

[![QAOA diagram](https://i.postimg.cc/sgryqr80/IMG-0202.jpg)](https://postimg.cc/XpQTBTnw)

Dependencies:

1) [Q# and Microsoft Quantum Libraries](https://docs.microsoft.com/en-us/quantum/language/)
2) [C# Accord Math library](http://accord-framework.net/docs/html/N_Accord_Math.htm)
