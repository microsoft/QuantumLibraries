#!/bin/env python
# -*- coding: utf-8 -*-
##
# packages.py: Abstraction to represent the list of packages.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

from distutils.version import LooseVersion
from typing import Iterable, Tuple

## LOGGING ##

import logging
logger = logging.getLogger(__name__)


## CLASSES ##

class Packages(object):
    """
    Represents the list of packages loaded into the current Q# session, and
    allows for adding new packages from NuGet.org or any other configured feeds.
    """

    def __init__(self, client):
        self._client = client

    def __iter__(self) -> Iterable[Tuple[str, LooseVersion]]:
        for pkg_spec in self._client.get_packages():
            name, version = pkg_spec.split("::", 1)
            yield name, LooseVersion(version)

    def __repr__(self) -> str:
        return repr(list(self))
    def __str__(self) -> str:
        return str(list(self))

    def add(self, package_name : str) -> None:
        """
        Adds a NuGet package with the given package name to the current Q#
        session, downloading the package from NuGet.org or any other configured
        feeds as necessary.
        """
        logger.info("Loading package: " + package_name)
        pkgs=self._client.add_package(package_name)
        logger.info("Loading complete: " + ';'.join(str(e) for e in pkgs))
