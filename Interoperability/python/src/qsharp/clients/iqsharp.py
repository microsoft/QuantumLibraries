#!/bin/env python
# -*- coding: utf-8 -*-
##
# iqsharp.py: Client for the IQ# Jupyter kernel.
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

class IQSharpError(RuntimeError):
    """
    Represents a Q# error passed by the IQ# kernel to the Python host.
    """
    def __init__(self, iqsharp_errors : List[str]):
        self.iqsharp_errors = iqsharp_errors
        error_msg = StringIO()
        error_msg.write("The Q# kernel raised the following errors:\n")
        error_msg.writelines([
            "    " + msg for msg in iqsharp_errors
        ])
        super().__init__(error_msg.getvalue())

class AlreadyExecutingError(IOError):
    """
    Raised when the IQ# client is already executing a command and cannot safely
    process an additional command.
    """
    pass

class IQSharpClient(object):
    kernel_manager = None
    kernel_client = None
    _busy : bool = False

    def __init__(self):
        self.kernel_manager = jupyter_client.KernelManager(kernel_name='iqsharp')

    ## Server Lifecycle ##

    def start(self):
        logger.info("Starting IQ# kernel...")
        self.kernel_manager.start_kernel(extra_arguments=["--user-agent", "qsharp.py"])
        self.kernel_client = self.kernel_manager.client()
        atexit.register(self.stop)

    def stop(self):
        logger.info("Stopping IQ# kernel...")
        try:
            self.kernel_manager.shutdown_kernel()
        except:
            pass

    def is_ready(self):
        try:
            result = self.component_versions(timeout=6)
            logger.info(f"Q# version\n{result}")
        except Exception as ex:
            logger.info('Exception while checking if IQ# is ready.', exc_info=ex)
            return
        return True

    def check_status(self):
        if not self.kernel_manager.is_alive():
            logger.debug("IQ# kernel is not running. Restarting.")
            self.start()

    ## Public Interface ##

    @property
    def busy(self) -> bool:
        return self._busy

    def compile(self, body):
        return self._execute(body)

    def get_available_operations(self) -> List[str]:
        return self._execute('%who', raise_on_stderr=False)

    def get_operation_metadata(self, name : str) -> Dict[str, Any]:
        return self._execute(f"?{name}")

    def get_workspace_operations(self) -> List[str]:
        return self._execute("%workspace")

    def reload(self) -> None:
        return self._execute(f"%workspace reload", raise_on_stderr=True)

    def add_package(self, name : str) -> None:
        return self._execute(f"%package {name}", raise_on_stderr=True)

    def get_packages(self) -> List[str]:
        return self._execute("%package", raise_on_stderr=False)

    def simulate(self, op, **params) -> Any:
        return self._execute_callable_magic('simulate', op, **params)

    def toffoli_simulate(self, op, **params) -> Any:
        return self._execute_callable_magic('toffoli', op, **params)

    def estimate(self, op, **params) -> Dict[str, int]:
        raw_counts = self._execute_callable_magic('estimate', op, **params)
        # Convert counts to ints, since they get turned to floats by JSON serialization.
        return {
            operation_name: int(count)
            for operation_name, count in raw_counts.items()
        }

    def component_versions(self, **kwargs) -> Dict[str, LooseVersion]:
        """
        Returns a dictionary from components of the IQ# kernel to their
        versions.
        """
        versions = {}
        def capture(msg):
            # We expect a display_data with the version table.
            if msg["msg_type"] == "display_data":
                data = unmap_tuples(json.loads(msg["content"]["data"]["application/json"]))
                for component, version in data["rows"]:
                    versions[component] = LooseVersion(version)
        self._execute("%version", output_hook=capture, **kwargs)
        return versions

    ## Internal-Use Methods ##

    def _execute_magic(self, magic : str, raise_on_stderr : bool = False, **params) -> Any:
        return self._execute(f'%{magic} {json.dumps(map_tuples(params))}', raise_on_stderr=raise_on_stderr)

    def _execute_callable_magic(self, magic : str, op, raise_on_stderr : bool = False, **params) -> Any:
        return self._execute_magic(f"{magic} {op._name}", raise_on_stderr=raise_on_stderr=, **params)

    def _execute(self, input, return_full_result=False, raise_on_stderr : bool = False, output_hook=None, **kwargs):
        logger.debug(f"sending:\n{input}")

        # make sure the server is still running:
        try:
            self.check_status()
        except:
            raise IQSharpError(["IQ# is not running."])

        results = []
        errors = []
        if output_hook is None:
            output_hook = self.kernel_client._output_hook_default
        def _output_hook(msg):
            if msg['msg_type'] == 'execute_result':
                results.append(msg)
            else:
                if raise_on_stderr and msg['msg_type'] == 'stream' and msg['content']['name'] == 'stderr':
                    errors.append(msg['content']['text'])
                else:
                    output_hook(msg)

        try:
            if self._busy:
                # Trying to execute while already executing can corrupt the
                # ordering of messages internally to ZeroMQ
                # (see https://github.com/Microsoft/QuantumLibraries/issues/69),
                # so we need to throw early rather than letting the problem
                # propagate to a Jupyter protocol error.
                raise AlreadyExecutingError("Cannot execute through the IQ# client while another execution is completing.")
            self._busy = True
            reply = self.kernel_client.execute_interactive(input, output_hook=_output_hook, **kwargs)
        finally:
            self._busy = False

        logger.debug(f"received:\n{reply}")

        # There should be either zero or one execute_result messages.
        if errors:
            raise IQSharpError(errors)
        if results:
            assert len(results) == 1
            content = results[0]['content']
            if 'application/json' in content['data']:
                obj = unmap_tuples(json.loads(content['data']['application/json']))
            else:
                obj = None
            return (obj, content) if return_full_result else obj
        else:
            return None
