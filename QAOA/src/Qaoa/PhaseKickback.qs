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
    /// ## phaseExponent
    /// Phase to be applied.
    /// ## controlQubits
    /// Qubits that are to aquire a phase.
    /// ## auxiliaryQubit
    /// auxiliary qubit.
    operation RunPhaseKickback(phaseExponent: Double, controlQubits: Qubit[]) : Unit is Adj + Ctl {
        using(auxiliaryQubit = Qubit()){
            within {
                ApplyToEachCA(
                    CNOT(_, auxiliaryQubit),
                    controlQubits
                );
            } apply {
                R(PauliZ, phaseExponent, auxiliaryQubit);
            }    
		}
	}
}