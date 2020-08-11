// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Measures all qubits and resets them.
    ///
    /// # Input
    /// ## qubits
    /// Qubits to be measured.
    ///
    /// # Output
    /// Results of the measurement.
   operation MeasureAllAndReset(qubits: Qubit[]) : Bool[] {
        let N = Length(qubits);
        mutable results = new Bool[N];
        for (i in 0..N-1) {
            set results w/= i <- (MResetZ(qubits[i]) == One);
        }
        return results;
    }

    /// # Summary
    /// Runs a phase kickback routine on an auxiliary qubit given indices of control qubits and a phase.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that are to aquire a phase.
    /// ## auxiliaryQubit
    /// auxiliary qubit.
    /// ## controlQubitsIndices
    /// List of indices of control qubits.
    /// ## phaseExponent
    /// Phase to be applied.
    operation RunPhaseKickback(qubits: Qubit[], auxiliaryQubit: Qubit, controlQubitsIndices: Int[], phaseExponent: Double) : Unit is Adj + Ctl {

        for(i in 0..Length(controlQubitsIndices)-1) {
            CNOT(qubits[controlQubitsIndices[i]], auxiliaryQubit);
		}

        R(PauliZ, phaseExponent, auxiliaryQubit);

        for(i in 0..Length(controlQubitsIndices)-1) {
            CNOT(qubits[controlQubitsIndices[i]], auxiliaryQubit);
		}
	}
}