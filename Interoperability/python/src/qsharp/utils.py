#!/bin/env python
# -*- coding: utf-8 -*-
##
# utils.py: Utilities internal to the qsharp package.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

import logging
logger = logging.getLogger(__name__)
from typing import Callable

## INTERNAL FUNCTIONS ##

def log_messages(data, action : Callable[[str], None] = logger.error):
    msgs = data['messages']
    for msg in msgs:
        action(msg)
