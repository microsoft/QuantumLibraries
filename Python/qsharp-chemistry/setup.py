#!/bin/env python
# -*- coding: utf-8 -*-
##
# setup.py: Installs Python host functionality for Q# chemistry library.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

## IMPORTS ##

import setuptools
import os

## VERSION INFORMATION ##
# Our build process sets the PYTHON_VERSION environment variable to a version
# string that is compatible with PEP 440, and so we inherit that version number
# here and propagate that to qsharp/version.py.
#
# To make sure that local builds still work without the same environment
# variables, we'll default to 0.0.0.1 as a development version.

version = os.environ.get('PYTHON_VERSION', '0.0.0.1')
is_conda = bool(os.environ.get('CONDA_BUILD', False))

with open('./qsharp_chemistry/version.py', 'w') as f:
    f.write(f'''# Auto-generated file, do not edit.
##
# version.py: Specifies the version of the qsharp-chemistry package.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##
__version__ = "{version}"
_user_agent_extra = "{"(conda)" if is_conda else ""}"
''')

## DESCRIPTION ##
# The long description metadata passed to setuptools is used to populate the
# PyPI page for this package. Thus, we'll generate the description by using the
# same README.md file that we use in the GitHub repo.

with open("./README.md", "r") as fh:
    long_description = fh.read()

## SETUPTOOLS INVOCATION ##

setuptools.setup(
    name="qsharp-chemistry",
    version=version,
    author="Microsoft",
    author_email="que-contacts@microsoft.com",
    description="Tools for interacting with the Q# chemistry library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/microsoft/QuantumLibraries",
    packages=setuptools.find_namespace_packages(include=["qsharp_chemistry.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'qsharp-core'
    ]
)
