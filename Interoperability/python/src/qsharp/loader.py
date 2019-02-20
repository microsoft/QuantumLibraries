#!/bin/env python
# -*- coding: utf-8 -*-
##
# loader.py: Support for exposing Q# namespaces as Python modules.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

import sys
from types import ModuleType
import importlib
from importlib.abc import MetaPathFinder, Loader
import qsharp

from typing import Optional, Any, Dict

import logging
logger = logging.getLogger(__name__)

class QSharpModuleFinder(MetaPathFinder):
    def find_module(self, full_name : str, path : Optional[str] = None) -> Loader:
        # We expose Q# namespaces as their own root-level packages.
        # E.g.:
        #     >>> import Microsoft.Quantum.Primitive as prim
        # Thus, we need to check if the full name is one that that we can
        # sensibly load before we proceed.
        ops = qsharp.get_available_operations_by_namespace()
        if full_name not in ops:
            # We may have been given part of the qualified name of a namespace.
            # E.g., if we try to import Microsoft.Quantum.Primitive, we'll
            # see calls with "Microsoft" and "Microsoft.Quantum" first.
            if not any(
                ns_name.startswith(full_name + ".")
                for ns_name in ops
            ):
                return None

        return QSharpModuleLoader()

class QSharpModuleLoader(Loader):
    def load_module(self, full_name : str):
        logger.debug(f"Trying to load {full_name} as a Q# namespace.")
        if full_name in sys.modules:
            return sys.modules[full_name]

        module = QSharpModule(full_name, full_name, self)

        # Register the new module.
        sys.modules.setdefault(full_name, module)
        return module

class QSharpCallable(object):
    _name : str
    def __init__(self, callable_name : str, source : str):
        self._name = callable_name
        self.source = source

    def __repr__(self) -> str:
        return f"<Q# callable {self._name}>"

    def simulate(self, **kwargs) -> Any:
        """
        Executes this function or operation on the QuantumSimulator target
        machine, returning its output as a Python object.
        """
        return qsharp.client.simulate(self, **kwargs)

    def estimate_resources(self, **kwargs) -> Dict[str, int]:
        return qsharp.client.estimate(self, **kwargs)

class QSharpModule(ModuleType):
    _qs_name : str

    def __init__(self, full_name : str, qs_name : str, loader : QSharpModuleLoader):
        super().__init__(full_name)
        self._qs_name = qs_name
        self.__file__ = f"qsharp:{qs_name}"
        self.__path__ = []
        self.__loader__ = loader

    def __getattr__(self, name):
        ops = qsharp.get_available_operations_by_namespace()
        if name in ops[self._qs_name]:
            return QSharpCallable(f"{self._qs_name}.{name}", "workspace")
        raise AttributeError(f"Q# namespace {self._qs_name} does not contain a callable {name}.")

    def __repr__(self) -> str:
        return f"<module '{self._qs_name}' (Q# namespace)>"
