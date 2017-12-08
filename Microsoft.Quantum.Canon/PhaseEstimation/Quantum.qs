// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Performs the quantum phase estimation algorithm for a given oracle U and targetState,
    /// reading the phase into a big-endian quantum register.
    ///
    /// # Input
    /// ## oracle
    /// An operation implementing U^m for given integer powers m.
    /// ## targetState
    /// A quantum register representing the state |φ〉 acted on by U. If |φ〉 is an
    /// eigenstate of U, U|φ〉 = e^{iφ} |φ〉 for φ ∈ [0, 2π) an unknown phase.
    /// ## controlRegister
    /// A big-endian representation integer register that can be used
    /// to control the provided oracle, and that will contain the a representation of φ following
    /// the application of this operation. The controlRegister is assumed to start in the initial
    /// state |00.0>, where the length of the register indicates the desired precision.
    operation QuantumPhaseEstimation(
		      oracle : DiscreteOracle, 
			  targetState : Qubit[],
			  controlRegister : BigEndian) : ()
    {
        body {
            let nQubits = Length(controlRegister);

			AssertAllZero(
				"`controlRegister` is expected to be in |0⟩⊗…⊗|0⟩",
				controlRegister, 1e-10);

            ApplyToEachCA(H, controlRegister);

            for (idxControlQubit in 0..(nQubits - 1)) {
                let control = controlRegister[idxControlQubit];
                let power = 2 ^ (nQubits - idxControlQubit - 1);
                (Controlled oracle)([control], (power, targetState));
            }

            QFT(controlRegister);
        }
		adjoint auto
		controlled auto
		controlled adjoint auto
    }
}
