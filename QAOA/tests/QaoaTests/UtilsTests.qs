// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.QAOA;

    @Test("QuantumSimulator")
    operation RunPhaseKickbackTest() : Unit {

        let numberOfQubits = 1;
        let auxiliaryQubits = 1;
        let controlQubitsIndices = [0];
        let phaseExponent = 0.5;

        using (qubits = Qubit[numberOfQubits+auxiliaryQubits]) {
            RunPhaseKickback(qubits[...1], qubits[1],controlQubitsIndices,phaseExponent);
            AssertQubit(Zero, Head(qubits));
            ResetAll(qubits);
        }
    }

     @Test("QuantumSimulator")
    operation RunPhaseKickbackOneControlQubitTest() : Unit {

        let numberOfQubits = 1;
        let auxiliaryQubits = 1;
        let controlQubitsIndices = [0];
        let phaseExponent = 0.5;
        
        let complexZero = Complex(0.685125,-0.174941);
        let complexOne = Complex(0.685125,0.174941);

        using (qubits = Qubit[numberOfQubits+auxiliaryQubits]) {
            H(qubits[0]);
            RunPhaseKickback([qubits[0]], qubits[1],controlQubitsIndices,phaseExponent);
            AssertQubitIsInStateWithinTolerance((complexZero,complexOne), qubits[0], 1E-05);
            ResetAll(qubits);
        }
    }

    @Test("QuantumSimulator")
    operation MeasureAllAndResetTest() : Unit {
        let numberOfQubits = 4;
        let complexOne = Complex(1.0, 0.0);
        let complexZero = Complex(0.0, 0.0);

        using (qubits = Qubit[numberOfQubits]) {
            X(qubits[0]);
            X(qubits[2]);
            let result = MeasureAllAndReset(qubits);

            Fact(result[0], "Expected |1> state.");
            Fact(not result[1], "Expected |0> state.");
            Fact(result[2], "Expected |1> state.");
            Fact(not result[3], "Expected |0> state.");

            AssertQubit(Zero, qubits[0]);
            AssertQubit(Zero, qubits[1]);
            AssertQubit(Zero, qubits[2]);
            AssertQubit(Zero, qubits[3]);

        }

    }

}
