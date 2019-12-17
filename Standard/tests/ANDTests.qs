// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;

    operation ANDTestHelper(polarity1 : Bool, polarity2 : Bool, gate : CCNOTop) : Unit {
        using ((control1, control2, target, output) = (Qubit(), Qubit(), Qubit(), Qubit())) {
            within {
                ApplyPauliFromBitString(PauliX, true, [polarity1, polarity2], [control1, control2]);
                gate::Apply(control1, control2, target);
            }
            apply {
                CNOT(target, output);
            }
            let expected = BoolAsResult(polarity1 and polarity2);
            if (MResetZ(output) != expected) {
                fail $"Expected output register to be {expected}";
            }
            if (IsResultOne(M(control1))) {
                fail "Expected first control register to be 0";
            }
            if (IsResultOne(M(control2))) {
                fail "Expected second control register to be 0";
            }
            if (IsResultOne(M(target))) {
                fail "Expected target register to be 0";
            }
        }
    }

    @Test("QuantumSimulator")
    operation ApplyANDTest() : Unit {
        for (p1 in [false, true]) {
            for (p2 in [false, true]) {
                for (op in [ApplyAND, ApplyANDLowDepth]) {
                    ANDTestHelper(p1, p2, CCNOTop(op));
                }
            }
        }
    }
}


