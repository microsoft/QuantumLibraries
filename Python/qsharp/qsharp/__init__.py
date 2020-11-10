# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.
import os
from distutils.util import strtobool

from qsharp.util import (
    init_qsharp,
    init_meta_path,
    compile,
    reload,
    get_available_operations,
    get_workspace_operations,
    get_available_operations_by_namespace,
    component_versions
)

from qsharp.clients.iqsharp import IQSharpError
from qsharp.types import Result, Pauli

# For debugging purposes only: set the environment variable 
# INIT_QSHARP=False to skip autoloading the qsharp client
if bool(strtobool(os.environ.get("INIT_QSHARP", "True"))):
    client, config, packages, projects = init_qsharp()
    init_meta_path()


__all__ = [
    'compile', 
    'reload',
    'get_available_operations', 
    'get_available_operations_by_namespace',
    'get_workspace_operations',
    'client',
    'config',
    'packages',
    'projects',
    'IQSharpError',
    'Result', 
    'Pauli'
]
