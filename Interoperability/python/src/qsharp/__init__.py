#!/bin/env python
# -*- coding: utf-8 -*-
##
# __init__.py: Root module for the qsharp package.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

"""

"""

## IMPORTS ##

import sys
from typing import List, Dict, Union
from collections import defaultdict

from qsharp.clients import _start_client
from qsharp.loader import QSharpCallable, QSharpModuleFinder
try:
    from qsharp.version import __version__
except:
    __version__ = "<unknown>"

## EXPORTS ##

__all__ = [
    'compile', 'get_available_operations', 'get_available_operations_by_namespace'
]

## FUNCTIONS ##

def compile(code : str) -> Union[QSharpCallable, List[QSharpCallable]]:
    """
    Given a string containing Q# source code, compiles it into the current
    workspace and returns one or more Q# callable objects that can be used to
    invoke the new code.

    :param code: A string containing Q# source code to be compiled.
    :returns: A list of callables compiled from `code`, or a callable if exactly
        one callable is found.
    """
    ops = [
        QSharpCallable(op, "snippets")
        for op in client.compile(code)
    ]
    if len(ops) == 1:
        return ops[0]

def get_available_operations() -> List[str]:
    """
    Returns a list containing the names of all operations and functions defined
    in the current workspace.
    """
    return client.get_available_operations()

def get_available_operations_by_namespace() -> Dict[str, List[str]]:
    """
    Returns a dictionary from namespaces to lists of operations and functions
    defined in the current workspace that are members of each namespace.
    """
    ops = get_available_operations()
    by_ns = defaultdict(list)

    for qualified_name in ops:
        idx_last_dot = qualified_name.rfind(".")
        ns_name = qualified_name[:idx_last_dot]
        op_name = qualified_name[idx_last_dot + 1:]

        by_ns[ns_name].append(op_name)

    return dict(by_ns.items())

## STARTUP ##

client = _start_client()

# Make sure that we're last on the meta_path so that actual modules are loaded
# first.
sys.meta_path.append(QSharpModuleFinder())
