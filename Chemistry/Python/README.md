# Q# Chemistry Library for Python #

The `qsharp-chemistry` package for Python provides interoperability with Microsoft Quantum Development Kit's Chemistry Lirary.

For details on how to get started with Python and Q#, please see the [Getting Started with Python guide](https://docs.microsoft.com/quantum/install-guide/python).

For details about the Quantum Chemistry Library, please see the [Introduction to the Quantum Chemistry Library article](https://docs.microsoft.com/en-us/quantum/user-guide/libraries/chemistry/) in our online documentation.

## Installing with Anaconda ##

If you use Anaconda or Miniconda, installing the `qsharp` package will automatically include all dependencies:

```bash
conda install -c quantum-engineering qsharp
```

## Installing from Source ##

If you'd like to contribute to or experiment with the Python interoperability feature, it may be useful to install from source rather than from the `qsharp` package on the Python Package Index (PyPI).
To do so, make sure that you are in the `Python` directory, and run `setup.py` with the `install` argument:

```bash
cd Chemistry/Python/
python setup.py install
```

## Building the `qsharp-chemistry` Package ##

The Python interoperability feature uses a standard `setuptools`-based packaging strategy.
To build a platform-independent wheel, run the setup script with `bdist_wheel` instead:

```bash
cd Chemistry/Python/
python setup.py bdist_wheel
```

By default, this will create a `qsharp-chemistry` wheel in `dist/` with the version number set to 0.0.0.1.
To provide a more useful version number, set the `PYTHON_VERSION` environment variable before running `setup.py`.
