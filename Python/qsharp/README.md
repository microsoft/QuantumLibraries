# Q# Interoperability for Python #

The `qsharp` package for Python provides interoperability with the Quantum Development Kit and with the Q# language, making it easy to simulate Q# operations and functions from within Python.

For details on how to get started with Python and Q#, please see the [Getting Started with Python guide](https://docs.microsoft.com/quantum/install-guide/python).

You can also try our [Quantum Computing Fundamentals](https://aka.ms/learnqc) learning path to get familiar with the basic concepts of quantum computing, build quantum programs, and identify the kind of problems that can be solved.

## Installing with Anaconda ##

If you use Anaconda or Miniconda, installing the `qsharp` package will automatically include all dependencies:

```bash
conda install -c quantum-engineering qsharp
```

## Installing from Source ##

If you'd like to contribute to or experiment with the Python interoperability feature, it may be useful to install from source rather than from the `qsharp` package on the Python Package Index (PyPI).

This package uses [namespace packages](https://www.python.org/dev/peps/pep-0382/) that are hosted in two metapackages: `qsharp-core` (part of [iqsharp](http://www.github.com/microsoft/iqsharp)) and `qsharp-chemistry` (part of [QuantumLibraries](http://www.github.com/microsoft/QuantumLibraries)). This means that to work in development mode, you need to install both.

To separately install `qsharp-chemistry` in development mode, run

```bash
cd Python/qsharp-chemistry
pip install -e .
```

If you also want to install `qsharp-core` in development mode, make sure to clone the [iqsharp](http://www.github.com/microsoft/iqsharp) repo and run

```bash
cd src/Python/qsharp-core
pip install -e .
```

To do so, make sure that you are in the `Python/qsharp` directory, and run `pip install -e .`.

```bash
cd Python/qsharp
pip install -e .
```

This will install `qsharp-core` and `qsharp-chemistry` from the PyPI if they were not previously installed.

## Building the `qsharp` Package ##

The Python interoperability feature uses a standard `setuptools`-based packaging strategy.
To build a platform-independent wheel, run the setup script with `bdist_wheel` instead:

```bash
cd Python/qsharp
python setup.py bdist_wheel
```

By default, this will create a `qsharp` wheel in `dist/` with the version number set to 0.0.0.1.
To provide a more useful version number, set the `PYTHON_VERSION` environment variable before running `setup.py`.

## Support and Q&A

If you have questions about the Quantum Development Kit and the Q# language, or if you encounter issues while using any of the components of the kit, you can reach out to the quantum team and the community of users in [Stack Overflow](https://stackoverflow.com/questions/tagged/q%23) and in [Quantum Computing Stack Exchange](https://quantumcomputing.stackexchange.com/questions/tagged/q%23) tagging your questions with **q#**.
