#!/bin/env python
# -*- coding: utf-8 -*-
##
# mock.py: Placeholder client useful in running unit tests and in hosted CI.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##


import subprocess
import time
import http.client
import atexit
import json
import sys
import urllib.parse
import os
import jupyter_client

from io import StringIO
from collections import defaultdict
from typing import List, Dict, Callable, Any
from pathlib import Path
from distutils.version import LooseVersion

from qsharp.utils import log_messages
from qsharp.serialization import map_tuples, unmap_tuples

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
        logger.info("Starting mock IQ# client...")
        atexit.register(self.stop)

    def stop(self):
        logger.info("Stopping mock IQ# client...")

    def is_ready(self):
        return True

    def check_status(self):
        pass

    ## Public Interface ##

    @property
    def busy(self) -> bool:
        return False

    def compile(self, body):
        return ["Workspace.Snippet.Example"]

    def get_available_operations(self) -> List[str]:
        return ["A.B.C", "A.B.D", "A.E.F"]

    def get_operation_metadata(self, name : str) -> Dict[str, Any]:
        return {}

    def get_workspace_operations(self) -> List[str]:
        return []

    def reload(self) -> None:
        return None

    def add_package(self, name : str) -> None:
        return self.packages.append(name)

    def get_packages(self) -> List[str]:
        return self.packages

    def simulate(self, op, **params) -> Any:
        return ()

    def estimate(self, op, **params) -> Dict[str, int]:
        return {
            "Depth": 13,
            "Width": 15
        }

    def component_versions(self, **kwargs) -> Dict[str, LooseVersion]:
        """
        Returns a dictionary from components of the IQ# kernel to their
        versions.
        """
        return {}
