// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Measurement {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Measures the given Pauli operator using an explicit scratch
    /// qubit to perform the measurement.
    ///
    /// # Input
    /// ## pauli
    /// A multi-qubit Pauli operator specified as an array of
    /// single-qubit Pauli operators.
    /// ## target
    /// Qubit register to be measured.
    ///
    /// # Output
    /// The result of measuring the given Pauli operator on
    /// the `target` register.
    operation MeasureWithScratch (pauli : Pauli[], target : Qubit[]) : Result {
        using (scratch = Qubit()) {
            H(scratch);

            for (idxPauli in IndexRange(pauli)) {
                let P = pauli[idxPauli];
                let src = target[idxPauli];

                if (P == PauliX) {
                    Controlled X([scratch], src);
                } elif (P == PauliY) {
                    Controlled Y([scratch], src);
                } elif (P == PauliZ) {
                    Controlled Z([scratch], src);
                }
            }

            H(scratch);
            return MResetZ(scratch);
        }
    }

    /// # Summary
    /// Given an array of multi-qubit Pauli operators, measures each using a specified measurement
    /// gadget, then returns the array of results.
    ///
    /// # Input
    /// ## paulis
    /// Array of multi-qubit Pauli operators to measure.
    /// ## target
    /// Register on which to measure the given operators.
    /// ## gadget
    /// Operation which performs the measurement of a given multi-qubit operator.
    ///
    /// # Output
    /// The array of results obtained from measuring each element of `paulis`
    /// on `target`.
    operation MeasurePaulis (paulis : Pauli[][], target : Qubit[], gadget : ((Pauli[], Qubit[]) => Result)) : Result[] {
        return ForEach(gadget(_, target), paulis);
    }

    /// # Summary
    /// Measures each qubit in a given array in the standard basis.
    /// # Input
    /// ## targets
    /// An array of qubits to be measured.
    /// # Output
    /// An array of measurement results.
    operation MultiM (targets : Qubit[]) : Result[] {
        return ForEach(M, targets);
    }

}
