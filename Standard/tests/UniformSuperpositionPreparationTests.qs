// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Preparation;


    @Test("QuantumSimulator")
    operation TestPrepareUniformSuperposition() : Unit {
        let nQubits = 5;
        use qubits = Qubit[nQubits];
        for nIndices in 1..2^nQubits {
            Message($"Testing nIndices {nIndices} on {nQubits} qubits");
            PrepareUniformSuperposition(nIndices, LittleEndian(qubits));

            let prob = 1.0 / IntAsDouble(nIndices);
            for stateIndex in 0..2^nQubits - 1 {
                AssertProbInt(stateIndex, stateIndex < nIndices ? prob | 0.0, LittleEndian(qubits), 1e-10);
            }

            ResetAll(qubits);
        }
    }

}
