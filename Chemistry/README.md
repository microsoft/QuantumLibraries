# Microsoft Quantum Chemistry Library #

This folder contains the C# and Q# sources used to implement the [Microsoft Quantum Chemistry library](https://docs.microsoft.com/azure/quantum/user-guide/libraries/chemistry).
Samples of how to use the library can be found in the Chemistry folder of the [Microsoft/Quantum repository](https://github.com/microsoft/Quantum/tree/main/Chemistry).

## Building and testing ##

The quantum chemistry library consists of two cross-platform project built using [.NET Core](https://docs.microsoft.com/en-us/dotnet/core/):

- [**DataModel.csproj**](https://github.com/microsoft/QuantumLibraries/tree/main/Chemistry/src/DataModel/DataModel.csproj): C# sources used to load, parse, and pre-compute Hamiltonians loaded from LIQ𝑈𝑖|〉 or Broombridge files.
- [**Runtime.csproj**](https://github.com/microsoft/QuantumLibraries/tree/main/Chemistry/src/Runtime/Runtime.csproj): Q# sources used to implement quantum chemistry simulation algorithms, given representations produced by the DataModel.

Once .NET Core is installed, you may build and run its tests by executing the following from a command line:

```bash
dotnet test tests
```

For more details about creating and running tests in Q#,
see the [Testing and debugging](https://docs.microsoft.com/azure/quantum/user-guide/testing-debugging)
section of the [developer's guide](https://docs.microsoft.com/azure/quantum).
