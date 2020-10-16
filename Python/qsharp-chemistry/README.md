# Q# Chemistry Library for Python #

The `qsharp-chemistry` package for Python provides interoperability with Microsoft Quantum Development Kit's Chemistry Library.

For details on how to get started with Python and Q#, please see the [Getting Started with Python guide](https://docs.microsoft.com/quantum/install-guide/python).

For details about the Quantum Chemistry Library, please see the [Introduction to the Quantum Chemistry Library article](https://docs.microsoft.com/quantum/user-guide/libraries/chemistry/) in our online documentation.

You can also try our [Quantum Computing Fundamentals](https://aka.ms/learnqc) learning path to get familiar with the basic concepts of quantum computing, build quantum programs, and identify the kind of problems that can be solved.

## Installing with Anaconda ##

If you use Anaconda or Miniconda, installing the `qsharp` package will automatically include all dependencies:

```bash
conda install -c quantum-engineering qsharp
```

## Installing from Source ##

If you'd like to contribute to or experiment with the Python interoperability feature, it may be useful to install from source rather than from the `qsharp-chemistry` package on the Python Package Index (PyPI).
To do so, make sure that you are in the `Python/qsharp-chemistry` directory, and run `setup.py` with the `install` argument:

```bash
cd Python/qsharp-chemistry
python setup.py install
```

## Building the `qsharp-chemistry` Package ##

The Python interoperability feature uses a standard `setuptools`-based packaging strategy.
To build a platform-independent wheel, run the setup script with `bdist_wheel` instead:

```bash
cd Python/qsharp-chemistry
python setup.py bdist_wheel
```

By default, this will create a `qsharp-chemistry` wheel in `dist/` with the version number set to 0.0.0.1.
To provide a more useful version number, set the `PYTHON_VERSION` environment variable before running `setup.py`.

## Support and Q&A

If you have questions about the Quantum Development Kit and the Q# language, or if you encounter issues while using any of the components of the kit, you can reach out to the quantum team and the community of users in [Stack Overflow](https://stackoverflow.com/questions/tagged/q%23) and in [Quantum Computing Stack Exchange](https://quantumcomputing.stackexchange.com/questions/tagged/q%23) tagging your questions with **q#**.
