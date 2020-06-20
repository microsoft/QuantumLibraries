// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;

    @Test("QuantumSimulator")
    operation CheckRepeatHIsNoOp() : Unit {
        AssertOperationsEqualReferenced(1,
            ApplyToEach(Repeat(H, 2, _), _),
            NoOp<Qubit[]>
        );

        AssertOperationsEqualReferenced(1,
            ApplyToEachC(RepeatC(H, 2, _), _),
            NoOp<Qubit[]>
        );

        AssertOperationsEqualReferenced(1,
            ApplyToEachA(RepeatA(H, 2, _), _),
            NoOp<Qubit[]>
        );

        AssertOperationsEqualReferenced(1,
            ApplyToEachCA(RepeatCA(H, 2, _), _),
            NoOp<Qubit[]>
        );
    }

}
