// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics as Diag;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    // Prepare single excitation
    operation PrepareSparseMultiConfigurationalState0Test () : Unit {
        let nQubits = 6;
        let intTest = [39];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5])];

        use qubits = Qubit[nQubits];
        PrepareSparseMultiConfigurationalState(NoOp, excitations, qubits);


        for idx in IndexRange(excitations) {
            AssertProbInt(intTest[idx], AbsD(1.0), LittleEndian(qubits), 1E-05);
        }

        ResetAll(qubits);
    }


    // Prepare multiple excitations with equal positive weights
    operation PrepareSparseMultiConfigurationalState1Test () : Unit {
        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5]), JordanWignerInputState((0.1, 0.0), [0, 2, 4]), JordanWignerInputState((0.1, 0.0), [3, 1])];

        use qubits = Qubit[nQubits];
        PrepareSparseMultiConfigurationalState(NoOp, excitations, qubits);

        for idx in IndexRange(excitations) {
            AssertProbInt(intTest[idx], expectedProb[idx], LittleEndian(qubits), 1E-05);
        }

        ResetAll(qubits);
    }


    // Prepare multiple excitations with unequal positive weights
    operation PrepareSparseMultiConfigurationalState2Test () : Unit {

        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [0.047619, 0.190476, 0.761905];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5]), JordanWignerInputState((0.2, 0.0), [0, 2, 4]), JordanWignerInputState((0.4, 0.0), [3, 1])];

        use qubits = Qubit[nQubits];
        PrepareSparseMultiConfigurationalState(NoOp, excitations, qubits);

        for idx in IndexRange(excitations) {
            AssertProbInt(intTest[idx], expectedProb[idx], LittleEndian(qubits), 1E-05);
        }

        ResetAll(qubits);
    }


    // Prepare multiple excitations with unequal complex weights
    operation PrepareSparseMultiConfigurationalState3Test () : Unit {
        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [0.047619, 0.190476, 0.761905];
        let p = [0.3, 0.9, 2.4];
        let excitations = [
            JordanWignerInputState(
                (0.1 * Cos(p[0]),
                0.1 * Sin(p[0])),
                [0, 1, 2, 5]
            ),

            JordanWignerInputState(
                (0.2 * Cos(p[1]),
                0.2 * Sin(p[1])),
                [0, 2, 4]
            ),

            JordanWignerInputState(
                (0.4 * Cos(p[2]), 0.4 * Sin(p[2])),
                [3, 1]
            )];

        use qubits = Qubit[nQubits];
        PrepareSparseMultiConfigurationalState(NoOp, excitations, qubits);

        for idx in IndexRange(excitations) {
            AssertProbInt(intTest[idx], expectedProb[idx], LittleEndian(qubits), 1E-05);
        }

        ResetAll(qubits);
    }


    // Prepare multiple excitations with complex weights
    operation PrepareSparseMultiConfigurationalState4Test () : Unit {
        let nQubits = 1;
        let intTest = [39, 21, 10];
        let phase = 2.453;
        let excitations = [JordanWignerInputState((0.1, 0.0), new Int[0]), JordanWignerInputState((0.1 * Cos(phase), 0.1 * Sin(phase)), [0])];

        use qubits = Qubit[nQubits];
        PrepareSparseMultiConfigurationalState(NoOp, excitations, qubits);
        Diag.AssertPhase(-phase / 2.0, qubits[0], 1E-09);
        ResetAll(qubits);
    }

}


