// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    operation CheckSetToBasisState(desired : Result) : Unit {
        using (q = Qubit()) {
            Ry(0.1234, q);
            SetToBasisState(desired, q);
            AssertQubit(desired, q);
            Reset(q);
        }
    }

    operation SetToBasisStateTest() : Unit {
        for (desired in [Zero, One]) {
            CheckSetToBasisState(desired);
        }
    }

}
