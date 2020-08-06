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
        let ancillaQubits = 1;
        let controlQubitsIndices = [0];
        let phaseExponent = 0.5;
        let complexZero = Complex(1.0, 0.0);
        let complexOne = Complex(0.0, 0.0);

        using (qubits = Qubit[numberOfQubits+ancillaQubits])
        {
            RunPhaseKickback(qubits[...1], qubits[1...],controlQubitsIndices,phaseExponent);
            AssertQubitIsInStateWithinTolerance((complexZero,complexOne), qubits[0], 1E-05);
            ResetAll(qubits);
        }
    }

     @Test("QuantumSimulator")
    operation RunPhaseKickbackTest2() : Unit {

        let numberOfQubits = 1;
        let ancillaQubits = 1;
        let controlQubitsIndices = [0];
        let phaseExponent = 0.5;
        
        let complexZero = Complex(0.685125,-0.174941);
        let complexOne = Complex(0.685125,0.174941);

        using (qubits = Qubit[numberOfQubits+ancillaQubits])
        {
            H(qubits[0]);
            RunPhaseKickback([qubits[0]], [qubits[1]],controlQubitsIndices,phaseExponent);
            AssertQubitIsInStateWithinTolerance((complexZero,complexOne), qubits[0], 1E-05);
            ResetAll(qubits);
        }
    }

    @Test("QuantumSimulator")
    operation MeasureAllAndResetTest() : Unit {
        let numberOfQubits = 4;
        let complexOne = Complex(1.0, 0.0);
        let complexZero = Complex(0.0, 0.0);

        using (qubits = Qubit[numberOfQubits])
        {
            X(qubits[0]);
            X(qubits[2]);
            let result = MeasureAllAndReset(qubits);

            EqualityFactB(result[0], true, "Expected |1> state.");
            EqualityFactB(result[1], false, "Expected |0> state.");
            EqualityFactB(result[2], true, "Expected |1> state.");
            EqualityFactB(result[3], false, "Expected |0> state.");

            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qubits[0], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qubits[1], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qubits[2], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qubits[3], 1E-06);

        }

    }

}
