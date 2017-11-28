---
title: "Measurement Example"
---

## Qubit Measurements ##

This example demonstrates the use of measurement operations to convert the state of one or more qubits into classical bits that can be used in the classical control code.

### Running the Example ###

Open the `QsharpLibraries.sln` solution in Visual Studio and set `Samples/Measurement/Measurement.csproj` as the startup project.
Press Start in Visual Studio to run the example.

### Manifest ###

- `Measurement/`
  - `Measurement.csproj`: Main C# project for the example.
  - `Measurement.qs`: Q# code preparing and measuring a few qubits.
  - `Driver.cs`: C# code to call the operations defined in Q#.
