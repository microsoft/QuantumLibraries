// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    function ZipTest() : () {
        let left = [1; 2; 101];
        let right = [PauliY; PauliI];
        let zipped = Zip(left, right);

        let (leftActual1, rightActual1) = zipped[0];
        if (leftActual1 != 1 || rightActual1 != PauliY) {
            fail $"Expected (1, PauliY), got ({leftActual1}, {rightActual1}).";
        }

        let (leftActual2, rightActual2) = zipped[1];
        if (leftActual2 != 2 || rightActual2 != PauliI) {
            fail $"Expected (2, PauliI), got ({leftActual2}, {rightActual2}).";
        }
    }

    function LookupTest() : () {
        let array = [1; 12; 71; 103];
        let fn = LookupFunction(array);
        AssertIntEqual(fn(0), 1, "fn(0) did not return array[0]");
        // Make sure we can call in random order!
        AssertIntEqual(fn(3), 103, "fn(3) did not return array[3]");
        AssertIntEqual(fn(2), 71, "fn(2) did not return array[2]");
        AssertIntEqual(fn(1), 12, "fn(1) did not return array[1]");
    }

    // TODO: Check Permute<'T>.
}
