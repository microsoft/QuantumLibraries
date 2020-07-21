// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Preparation;

    /// # Summary
    /// Given an operation, displays diagnostics about
    /// the operation that are made available by the current
    /// execution target.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits on which the given operation acts.
    /// ## op
    /// The operation that is to be diagnosed.
    ///
    /// # Remarks
    /// Calling this operation has no observable effect from within
    /// Q#. The exact diagnostics that are displayed, if any, are
    /// dependent on the current execution target and editor environment.
    /// For example, when used on the full-state quantum simulator,
    /// a unitary matrix used to represent `op` is displayed.
    operation DumpOperation(nQubits : Int, op : (Qubit[] => Unit is Adj))
    : Unit {
        using ((reference, target) = (Qubit[nQubits], Qubit[nQubits])) {
            // The operation provided could be a partial application of
            // another operation, such that there could be an observable
            // effect of dumping this operation unless we undo preparing the
            // Choi state.
            within {
                PrepareChoiStateA(op, reference, target);
            } apply {
                DumpReferenceAndTarget(reference, target);
            }
        }
    }

}
