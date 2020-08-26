// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;

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

        within {
            ApplyToEachCA(
                CNOT(_, auxiliaryQubit),
                Subarray(controlQubitsIndices, qubits)
            );
        } apply {
            R(PauliZ, phaseExponent, auxiliaryQubit);
        }
	}
}