// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Diagnostics;

    function NotTest() : Unit {
        Fact(Not(false), "Not returned wrong output.");
    }

    function AndTest() : Unit {
        Fact(not And(false, false), "And returned wrong output.");
        Fact(not And(false, true), "And returned wrong output.");
        Fact(not And(true, false), "And returned wrong output.");
        Fact(And(true, true), "And returned wrong output.");
    }

    function OrTest() : Unit {
        Fact(not Or(false, false), "Or returned wrong output.");
        Fact(Or(false, true), "Or returned wrong output.");
        Fact(Or(true, false), "Or returned wrong output.");
        Fact(Or(true, true), "Or returned wrong output.");
    }

    function XorTest() : Unit {
        Fact(not Xor(false, false), "Xor returned wrong output.");
        Fact(Xor(false, true), "Xor returned wrong output.");
        Fact(Xor(true, false), "Xor returned wrong output.");
        Fact(not Xor(true, true), "Xor returned wrong output.");
    }

    function ConditionedTest() : Unit {
        EqualityFactI(Conditioned(true, 42, -1), 42, "Conditioned returned wrong output.");
        EqualityFactL(Conditioned(false, 42L, -1L), -1L, "Conditioned returned wrong output.");
    }

}
