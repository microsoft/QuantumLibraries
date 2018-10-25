# Quantum Chemistry Library Tests #

Unit tests for the Microsoft Quantum Chemistry library. All files in this library are intended to eventually be open-sourced under the MIT license, like the other components of github.com/microsoft/quantum.

## Structure of the tests

- **[DataModelTests](DataModelTests/)**:
    Contains tests for the C# component of handling the chemistry Hamiltonian and loading them from files. 

- **[ChemistryTests](ChemistryTests/)**:
    Contains tests for isolated Q# component of simulating chemistry Hamiltonians.

- **[SystemTests](SystemTests/)**:
    Contains tests for simulating chemistry Hamiltonians that are output by the DataModel library. 

- **[CanonAdditionsTests](CanonAdditionsTests/)**:
    These tests of the Q# canon additions focus currently mainly on LCU techniques, in particular on qubitization. 
    