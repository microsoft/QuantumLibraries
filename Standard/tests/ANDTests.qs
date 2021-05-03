// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.ANDTests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Measurement;

    operation AndTestHelper(polarity1 : Bool, polarity2 : Bool, gate : CCNOTop) : Unit {
        use control1 = Qubit();
        use control2 = Qubit();
        use target = Qubit();
        use output = Qubit();
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

    operation ControlledAndTestHelper(polarities : Bool[], gate : ((Qubit, Qubit, Qubit) => Unit is Adj+Ctl)) : Unit {
        let numControls = Length(polarities);
        use controls = Qubit[numControls];
        use target = Qubit();
        use output = Qubit();
        within {
            ApplyPauliFromBitString(PauliX, true, polarities, controls);
            Controlled gate(controls[0..numControls - 3], (controls[numControls - 2], controls[numControls - 1], target));
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

    @Test("QuantumSimulator")
    @Test("ToffoliSimulator")
    operation ApplyAndTest() : Unit {
        for p1 in [false, true] {
            for p2 in [false, true] {
                for op in [ApplyAnd, ApplyLowDepthAnd] {
                    AndTestHelper(p1, p2, CCNOTop(op));
                }
            }
        }
    }

    @Test("QuantumSimulator")
    @Test("ToffoliSimulator")
    operation ControlledApplyAndTest() : Unit {
        for numControls in 3..5 {
            for assignment in 0..2^numControls - 1 {
                ControlledAndTestHelper(IntAsBoolArray(assignment, numControls), ApplyAnd);
            }
        }
        for numControls in 3..4 {
            for assignment in 0..2^numControls - 1 {
                ControlledAndTestHelper(IntAsBoolArray(assignment, numControls), ApplyLowDepthAnd);
            }
        }
    }
}
