# Quantum Chemistry Library Tests #

Unit tests for the [Microsoft Quantum Chemistry library](https://docs.microsoft.com/en-us/quantum/libraries/chemistry/).

## Structure of the tests

- **[DataModelTests](DataModelTests/)**:
    Contains tests for the C# component of handling the chemistry Hamiltonian and loading them from files.

- **[ChemistryTests](ChemistryTests/)**:
    Contains tests for isolated Q# component of simulating chemistry Hamiltonians.

- **[SystemTests](SystemTests/)**:
    Contains tests for simulating chemistry Hamiltonians that are output by the DataModel library.
