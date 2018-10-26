# Microsoft Quantum Chemistry Library. 

C# and Q# sources used to implement the Microsoft Quantum Chemistry library. Samples and tests included.

## Outline of features
- Load a chemistry Hamiltonian from file. Supported format are:
    - LiQui|> schema
    - YAML schema
- Pass loaded Hamiltonian to a variety of simulation algorithm.
    - Trotter simulation
    - Optimized Trotter simulation
    - Qubitization with minimal qubit overhead
    - Qubitization with minimal T-gate overhead
- Samples for learning how to use the library, and also for performing quantum simulations of ground state energy estimation and obtaining resource estimates of simulation algorithms.

## Verify Installation
- If using Microsoft Visual Studio:
1. Open 'Microsoft.Quantum.Chemistry.sln'.
2. Select Samples/1 - MolecularHydrogen/MolecularHydrogenGUI as the StartUp project.
3. Press F5 to run the molecular Hydrogen quantum phase estimation demo.

- If using Windows command line:
1. Go to [RunSimulation](../Samples/Chemistry/MolecularHydrogenGUI).
2. Enter 'dotnet run' to run molecular Hydrogen quantum phase estimation demo.


## Structure of the Library

- **[DataModel](src/DataModel/)**:
    This is a C# library that handles loading a chemistry Hamiltonian from file into a standard format. This library also handles classical preprocessing & optimization of the standard format to a format specialized for consumption by various Q# simulation algorithms.

- **[Chemistry](src/Chemistry/)**:
    This Q# library contains methods specific to simulating chemistry Hamiltonians that are output by the DataModel library.

- **[CanonAdditions](src/CanonAdditions/)**:
    This Q# library contains methods that will be eventually merged with the Microsoft.Quantum.Canon library.

## Structure of the tests

- **[DataModelTests](tests/DataModelTests/)**:
    Contains tests for the C# component of handling the chemistry Hamiltonian and loading them from files. 

- **[ChemistryTests](tests/ChemistryTests/)**:
    Contains tests for isolated Q# component of simulating chemistry Hamiltonians.

- **[SystemTests](tests/SystemTests/)**:
    Contains tests for simulating chemistry Hamiltonians that are output by the DataModel library. 

## Structure of the samples
See the README in [Chemistry library samples](../Samples/Chemistry).