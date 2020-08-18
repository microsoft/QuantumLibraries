// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {

    /// # Summary
    /// Between a call to this operation and its adjoint, asserts that
    /// a given operation is called at most a certain number of times.
    ///
    /// # Input
    /// ## nTimes
    /// The maximum number of times that `op` may be called.
    /// ## op
    /// An operation whose calls are to be restricted.
    /// ## message
    /// A message to be displayed upon failure.
    ///
    /// # Example
    /// The following snippet will fail when executed on machines which
    /// support this diagnostic:
    /// ```Q#
    /// using (register = Qubit[4]) {
    ///     within {
    ///         AllowAtMostNCallsCA(3, H, "Too many calls to H.");
    ///     } apply {
    ///         // Fails since this calls H four times, rather than the
    ///         // allowed maximum of three.
    ///         ApplyToEach(H, register);
    ///     }
    /// }
    /// ```
    ///
    /// # Remarks
    /// This operation may be replaced by a no-op on targets which do not
    /// support it.
    operation AllowAtMostNCallsCA<'TInput, 'TOutput>(
        nTimes : Int, op : ('TInput => 'TOutput is Adj + Ctl),
        message : String
    )
    : Unit is Adj {
    }

    /// # Summary
    /// Between a call to this operation and its adjoint, asserts that
    /// at most a given number of additional qubits are allocated with
    /// using statements.
    ///
    /// # Input
    /// ## nQubits
    /// The maximum number of qubits that may be allocated.
    /// ## message
    /// A message to be displayed upon failure.
    ///
    /// # Example
    /// The following snippet will fail when executed on machines which
    /// support this diagnostic:
    /// ```Q#
    /// within {
    ///     AllowAtMostNQubits(3, "Too many qubits allocated.");
    /// } apply {
    ///     // Fails since this allocates four qubits.
    ///     using (register = Qubit[4]) {
    ///     }
    /// }
    /// ```
    ///
    /// # Remarks
    /// This operation may be replaced by a no-op on targets which do not
    /// support it.
    operation AllowAtMostNQubits(nQubits : Int, message : String) : Unit is Adj {
    }

}
