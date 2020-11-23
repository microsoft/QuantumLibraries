# Microsoft Quantum QAOA Library

This library provides a hybrid quantum-classical algorithm for solving optimization problems. 
It includes a Q# implementation of the Quantum Approximate Optimization Algorithm ([QAOA](https://arxiv.org/abs/1411.4028)) together with a classical optimizer in C#.
The classical optimizer uses a quantum objective function to choose hopefully optimal parameters for running the QAOA.
The quantum objective function is currently evaluated on a simulator backend provided by Q#.

## Building and testing ##

The quantum QAOA library consists of a cross-platform project built using [.NET Core](https://docs.microsoft.com/en-us/dotnet/core/):

- **Qaoa.csproj**: Q# and C# sources used to encode combinatorial optimization problems and run QAOA and hybrid QAOA algorithms.

The high-level diagram of the implementation (notation comes from the [QAOA paper](https://arxiv.org/abs/1411.4028)):

[![QAOA diagram](https://i.postimg.cc/sgryqr80/IMG-0202.jpg)](https://postimg.cc/XpQTBTnw)

Once .NET Core is installed, you may build and run its tests by executing the following from a command line:

```bash
dotnet test tests
```

For more details about creating and running tests in Q#,
see the [Testing and debugging](https://docs.microsoft.com/quantum/quantum-techniques-testinganddebugging)
section of the [developer's guide](https://docs.microsoft.com/quantum).

Dependencies:

1) [Q# and Microsoft Quantum Libraries](https://docs.microsoft.com/en-us/quantum/language/)
2) [C# Accord Math library](http://accord-framework.net/docs/html/N_Accord_Math.htm)
