---
title: "Ising Trotter Sample"
---

## Ising Trotter Sample ##

This sample walks through constructing the time-evolution operator for the Ising model using the Trotterization library feature. This time-evolution operator is applied to investigate spin relaxation.

### Running the Sample ###

Open the `QsharpLibraries.sln` solution in Visual Studio and set the .csproj file in the manifest as the startup project.
Press Start in Visual Studio to run the sample.

### Manifest ###

- `IsingTrotterSample/`
  - `IsingTrotterSample.csproj`: Main C# project for the sample.
  - `IsingTrotter.qs`: Q# code implementing quantum operations for this sample.
  - `Program.cs`: C# code to interact with and print out results of the Q# operations for this sample.
