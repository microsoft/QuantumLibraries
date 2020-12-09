#!/bin/env python
# -*- coding: utf-8 -*-
##
# setup.py: Installs Python host functionality for Q# chemistry library.
##
# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
##

## IMPORTS ##

import setuptools
import os

version = os.environ.get('PYTHON_VERSION', '0.0.0.1')
is_conda = bool(os.environ.get('CONDA_BUILD', False))

## DESCRIPTION ##
# Chemistry tools for the Microsoft Quantum Development Kit

with open("./README.md", "r") as fh:
    long_description = fh.read()

## SETUPTOOLS INVOCATION ##

setuptools.setup(
    name="qdk_chemistry",
    version=version,
    author="Microsoft",
    author_email="que-contacts@microsoft.com",
    description="Chemistry library for the QDK.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/microsoft/QuantumLibraries",
    packages=setuptools.find_namespace_packages(include=["qdk*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'qsharp',
        'jupyter_jsmol',
        'networkx',
        'varname'
    ]
)
