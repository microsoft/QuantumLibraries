// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Measurement {
    open Microsoft.Quantum.Intrinsic;
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

    /// # Summary
    /// Measures a register of qubits and returns true if it is in the all-zeros state or false otherwise.
    ///
    /// # Description
    /// Given a register of qubits, checks if that register is in the state
    /// $\ket{00 \cdots 0}$ by performing a computational basis (i.e.:
    /// `PauliZ`) measurement on each individual qubit in the register.
    ///
    /// # Input
    /// ## register
    /// The register of qubits to be measured.
    ///
    /// # Output
    /// `true` if and only if the register is measured to be in the
    /// $\ket{00 \cdots 0}$ state.
    ///
    /// # Remarks
    /// This operation does not reset its qubits, but projects them to a
    /// computational basis state.
    operation MeasureIfAllZeros(register : Qubit[]) : Bool {
        return All(IsResultZero, ForEach(M, register));
    }

}
