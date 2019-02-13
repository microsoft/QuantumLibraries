# Python Interoperability for Q# #

The `qsharp` package for Python provides interoperability with the Quantum Development Kit and with the Q# language, making it easy to simulate Q# operations and functions from within Python.

## Installation ##

To get started with the `qsharp` package, you'll need a couple prerequisites:

- .NET Core SDK 2.1.3 or later,
- Python 3.6 or later,
- Jupyter, and
- The IQ# kernel for Jupyter.

To install the .NET Core SDK, please follow the [Hello World tutorial](https://dotnet.microsoft.com/learn/dotnet/hello-world-tutorial/intro) provided with .NET Core.

To install Python and Jupyter, we recommend using the Anaconda distribution of Python.
Please see https://www.anaconda.com/distribution/ for more details.
On other distributions of Python, the Jupyter platform can be installed using `pip install jupyter`.

Finally, the IQ# kernel for Jupyter can be installed using the `dotnet` command line tool:

```
dotnet tool install -g Microsoft.Quantum.IQSharp
dotnet iqsharp install
```

For more information, please see the complete [install guide](https://docs.microsoft.com/quantum/install-guide/).

# Usage

Create one or more files with a `.qs` extension with the quantum operations you want to execute.
The `qsharp` package automatically detects and tries to compile all the files under the current directory that have the  `.qs` extension.

To call a Q# operation from Python, first import `qsharp`:
```python
import qsharp
```

After this, Q# namespaces can be imported as Python packages, for example:
```python
from Microsoft.Quantum.Python import HelloQ, HelloAgain
```

Once imported, to simulate a Q# operation invoke it's `simulate` method:
```python
HelloQ.simulate()
```

If the Q# operation expects parameters, include them as named parameters to the `simulate` method:
```python
HelloAgain.simulate(count=3, name="Alice")
```

If the Q# operation returns a value, the corresponding value is returned from the `simulate` method.
```python
r = HelloAgain.simulate(count=3, name="Alice")
print("HelloAgain result: ", r)
```

On top of simulation, you can also do quantum resources estimation including 
the count of primitive operations used by the algorithm and the number of required qubits.
For this, invoke the `trace` method on the operation:
```python
r = HelloAgain.trace(count=5, name="Counting")
```

You may use the `qsharp.print_tracer_counts` method to print the trace results to the console:
```python
qsharp.print_tracer_counts(r)
```


You can now also compile Q# operations on the fly from Python
and simulate them.
To create an operation on the fly, the `qsharp` module exports a `compile` method which receives two parmaters
* id: a unique identifier of this snippet.
* source: a valid Q# code snippet with the operation definition, for example:
```python
hello = qsharp.compile("snippet_1", """
    operation HelloQ() : Result
    {
        Message($"Hello from quantum world!"); 
        return Zero;
    }
""")
```

If successful, `compile` returns a Q# operation that can now be simulated or traced:
```
r = hello.simulate()
```

You may call `compile` multiple times, and may refer to operations previously defined in other snippets. 
Calling `compile` using a previous snippet_id updates the corresponding definition.
