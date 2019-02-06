#!/bin/env python
# -*- coding: utf-8 -*-
##
# driver.py: Example of a script that calls into Q# from Python.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

# To start, simply import the qsharp module.
import qsharp

print(qsharp.get_available_operations())

# After importing the qsharp module, Q# namespaces can be imported like any
# other Python packages.
from Microsoft.Prototypes.Python import HelloQ, HelloAgain, HelloTuple

# All Q# operations defined in any .qs file inside the same working directory
# is automatically identified.
# You can simulate any of them simply by calling the `simulate` method
# on Q# functions and operations once they've been imported.
r = HelloQ.simulate()
print("HelloQ result: ", r)
print("")

# If the operation receives parameters, just include them as named parameters
# of the same simulate method.
r = HelloAgain.simulate(count=3, name="Counting")
print("HelloAgain result: ", r)
print("")

# All built-in types are currently supported:
r = HelloTuple.simulate(count=2, tuples=[(1, "uno"), (2, "dos"), (3, "tres"), (4, "cuatro")])
print("HelloTuple result: ", r)
print("")

# On top of simulation, you can also do quantum resources estimation including
# the count of primitive operations used by the algorithm and the number of required qubits.
# For this, invoke the `estimate` method on the operation:
r = HelloAgain.estimate(count=5, name="Counting")
print()

# # You may use the `print_tracer_counts` function `qsharp` to print the results to the console:
qsharp.print_tracer_counts(r)


# You can now also compile Q# operations on the fly from Python
# and simulate them.
# To create an operation on the fly, call the compile method which receives two parmaters
#  * id: a unique identifier of this snippet.
#  * source: a valid Q# code snippet with the operation definition, for example:
hello = qsharp.compile("""
    operation HelloQ() : Result
    {
        Message($"Hello from quantum world!"); 
        return Zero;
    }
""")

# if successful, `compile` returns a Q# operation that can now be simulated or traced:
r = hello.simulate()
print("First snippet: ", r)
print()

# You may call `compile` multiple times, and may refer to operations previously defined in other snippets:
call_hello = qsharp.compile("""
    operation CallHello() : Bool {
        Message("Calling HelloQ");
        if (HelloQ() == One) {
            return true;
        } else {
            return false;
        }
    }
""")
r = call_hello.simulate()
print("Second snippet: ", r)
print()

# Calling `compile` using a previous snippet id updates the corresponding definition
hello = qsharp.compile("""
    operation HelloQ() : Result
    {
        Message("MSG1");
        Message("MSG2");

        return One;
    }
""")
r = call_hello.simulate()
print("Revised snippet: ", r)
print()
