// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Measurement;

    operation AndTestHelper(polarity1 : Bool, polarity2 : Bool, gate : CCNOTop) : Unit {
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
            AssertAllZero([control1, control2, target]);
        }
    }

    operation ControlledAndTestHelper(polarities : Bool[]) : Unit {
        let numControls = Length(polarities);
        using ((controls, target, output) = (Qubit[numControls], Qubit(), Qubit())) {
            within {
                ApplyPauliFromBitString(PauliX, true, polarities, controls);
                Controlled ApplyAnd(controls[0..numControls - 3], (controls[numControls - 2], controls[numControls - 1], target));
            }
            apply {
                CNOT(target, output);
            }
            let expected = BoolAsResult(All(EqualB(true, _), polarities));
            if (MResetZ(output) != expected) {
                fail $"Expected output register to be {expected}";
            }
            AssertAllZero(controls + [target]);
        }
    }

    @Test("QuantumSimulator")
    operation ApplyAndTest() : Unit {
        for (p1 in [false, true]) {
            for (p2 in [false, true]) {
                for (op in [ApplyAnd, ApplyLowDepthAnd]) {
                    AndTestHelper(p1, p2, CCNOTop(op));
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation ControlledApplyAndTest() : Unit {
        for (numControls in 3..5) {
            for (assignment in 0..2^numControls - 1) {
                ControlledAndTestHelper(IntAsBoolArray(assignment, numControls));
            }
        }
    }
}


