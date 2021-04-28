// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    operation CheckSetToBasisStateHelper(desired : Result) : Unit {
        use q = Qubit();
        Ry(0.1234, q);
        SetToBasisState(desired, q);
        AssertQubit(desired, q);
        Reset(q);
    }

    @Test("QuantumSimulator")
    operation CheckSetToBasisState() : Unit {
        for desired in [Zero, One] {
            CheckSetToBasisStateHelper(desired);
        }
    }

}
