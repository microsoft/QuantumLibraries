#!/bin/env python
# -*- coding: utf-8 -*-
##
# loader.py: Support for exposing Q# namespaces as Python modules.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

import sys
from types import ModuleType, new_class
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
        #     >>> import Microsoft.Quantum.Intrinsic as mqi
        # Thus, we need to check if the full name is one that that we can
        # sensibly load before we proceed.

        # To check the full name, we ask the client rather than going through
        # the public API for the qsharp package, so that we can check if the
        # client is currently busy. This can happen if anything below us in
        # meta_path needs to handle an import during an execute; this is the
        # case when ZeroMQ needs to import additional functionality from a
        # Cython module to handle a message.
        # See https://github.com/Microsoft/QuantumLibraries/issues/69 for an
        # example of this failure modality.

        # If the client is busy, we'll want to forego this request to find a
        # module and return None early.
        if qsharp.client.busy:
            return None

        # At this point, we should be safe to rely on the public API again.
        ops = qsharp.get_available_operations_by_namespace()

        if full_name not in ops:
            # We may have been given part of the qualified name of a namespace.
            # E.g., if we try to import Microsoft.Quantum.Intrinsic, we'll
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
            op_cls = new_class(name, (QSharpCallable, ))

            # Copy over metadata from the operation's header.
            metadata = qsharp.client.get_operation_metadata(f"{self._qs_name}.{name}")
            op_cls.__doc__ = metadata.get('documentation', '')
            op_cls.__file__ = metadata.get('source', None)
            return op_cls(f"{self._qs_name}.{name}", "workspace")
        raise AttributeError(f"Q# namespace {self._qs_name} does not contain a callable {name}.")

    def __repr__(self) -> str:
        return f"<module '{self._qs_name}' (Q# namespace)>"
