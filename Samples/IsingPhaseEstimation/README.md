---
title: "Ising Phase Estimation Sample"
---

## Ising Phase Estimation Sample ##

This sample adiabatically prepares the ground state of the Ising model Hamiltonian, and then perform phase estimation to obtain an estimate of the ground state energy.

### Running the Sample ###

Open the `QsharpLibraries.sln` solution in Visual Studio and set the .csproj file in the manifest as the startup project.
Press Start in Visual Studio to run the sample.

### Manifest ###

- `IsingPhaseEstimation/`
  - `IsingPhaseEstimationSample.csproj`: Main C# project for the sample.
  - `IsingPhaseEstimation.qs`: Q# code implementing quantum operations for this sample.
  - `Program.cs`: C# code to interact with and print out results of the Q# operations for this sample.

## Note ##

This sample builds on results in `AdiabaticIsingSample` and uses techniques introduced in `PhaseEstimationSample`.
We suggest reading those samples before continuing.

