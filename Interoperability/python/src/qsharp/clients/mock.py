#!/bin/env python
# -*- coding: utf-8 -*-
##
# mock.py: Placeholder client useful in running unit tests and in hosted CI.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

import atexit
import json

from typing import List, Dict, Callable, Any
from distutils.version import LooseVersion

## LOGGING ##

import logging
logger = logging.getLogger(__name__)

## CLASSES ##

class MockClient(object):
    packages: List[str]

    def __init__(self):
        self.packages = []

    ## Server Lifecycle ##

    def start(self):
        logger.debug("MockClient.start called.")
        atexit.register(self.stop)

    def stop(self):
        logger.debug("MockClient.stop called.")

    def is_ready(self):
        logger.debug("MockClient.is_ready called.")
        return True

    def check_status(self):
        logger.debug("MockClient.check_status called.")

    ## Public Interface ##

    @property
    def busy(self) -> bool:
        logger.debug("MockClient.busy accessed.")
        return False

    def compile(self, body):
        logger.debug(f"MockClient.compile called with body:\n{body}")
        return ["Workspace.Snippet.Example"]

    def get_available_operations(self) -> List[str]:
        logger.debug("MockClient.get_available_operations called.")
        return ["A.B.C", "A.B.D", "A.E.F"]

    def get_operation_metadata(self, name : str) -> Dict[str, Any]:
        logger.debug(f"MockClient.get_operation_metadata called with name {name}.")
        return {}

    def get_workspace_operations(self) -> List[str]:
        logger.debug("MockClient.get_workspace_operations called.")
        return []

    def reload(self) -> None:
        logger.debug("MockClient.reload called.")
        return None

    def add_package(self, name : str) -> None:
        logger.debug(f"MockClient.add_package called with name {name}.")
        return self.packages.append(name)

    def get_packages(self) -> List[str]:
        logger.debug("MockClient.get_packages called.")
        return self.packages

    def simulate(self, op, **params) -> Any:
        logger.debug(f"MockClient.simulate called with operation {op} and params:\n{params}")
        return ()

    def estimate(self, op, **params) -> Dict[str, int]:
        logger.debug(f"MockClient.estimate called with operation {op} and params:\n{params}")
        return {
            "Depth": 13,
            "Width": 15
        }

    def component_versions(self, **kwargs) -> Dict[str, LooseVersion]:
        """
        Returns a dictionary from components of the IQ# kernel to their
        versions.
        """
        logger.debug(f"MockClient.component_versions called with keyword arguments:\n{kwargs}")
        return {}
