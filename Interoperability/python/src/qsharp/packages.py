#!/bin/env python
# -*- coding: utf-8 -*-
##
# packages.py: Abstraction to represent the list of packages.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

from distutils.version import LooseVersion

## CLASSES ##

class Packages(object):
    def __init__(self, client):
        self._client = client

    def __iter__(self):
        for pkg_spec in self._client.get_packages():
            name, version = pkg_spec.split("::", 1)
            yield name, LooseVersion(version)

    def __repr__(self):
        return repr(list(self))
    def __str__(self):
        return str(list(self))

    def add(self, package_name : str) -> None:
        self._client.add_package(package_name)
