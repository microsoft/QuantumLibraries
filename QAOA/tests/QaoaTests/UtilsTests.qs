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
        let complexOne = Complex(1.0, 0.0);
        let complexZero = Complex(0.0, 0.0);

        using (qubits = Qubit[numberOfQubits+ancillaQubits])
        {
            RunPhaseKickback(qubits[...1], qubits[1...],controlQubitsIndices,phaseExponent);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qubits[0], 1E-06);
            ResetAll(qubits);
        }
    }

    @Test("QuantumSimulator")
    operation MeasureAllAndResetTest() : Unit {


    }

}
