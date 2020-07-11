# Python Interoperability for Q# #

The `qsharp` package for Python provides interoperability with the Quantum Development Kit and with the Q# language, making it easy to simulate Q# operations and functions from within Python.

For details on how to get started with Python and Q#, please see the [Getting Started with Python guide](https://docs.microsoft.com/quantum/install-guide/python).

## Installing with Anaconda ##

If you use Anaconda or Miniconda, installing the `qsharp` package will automatically include all dependencies:

```bash
conda install -c quantum-engineering qsharp
```

## Installing from Source ##

If you'd like to contribute to or experiment with the Python interoperability feature, it may be useful to install from source rather than from the `qsharp` package on the Python Package Index (PyPI).
To do so, make sure that you are in the `src` directory, and run `setup.py` with the `install` argument:

```bash
cd iqsharp/src/Python/
python setup.py install
```

## Building the `qsharp` Package ##

The Python interoperability feature uses a standard `setuptools`-based packaging strategy.
To build a platform-independent wheel, run the setup script with `bdist_wheel` instead:

```bash
cd iqsharp/src/Python/
python setup.py bdist_wheel
```

By default, this will create a `qsharp` wheel in `dist/` with the version number set to 0.0.0.1.
To provide a more useful version number, set the `ASSEMBLY_VERSION` environment variable before running `setup.py`.
