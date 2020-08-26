// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Qaoa;

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
}
