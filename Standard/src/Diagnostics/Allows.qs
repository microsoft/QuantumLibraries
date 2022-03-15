// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {

    /// # Summary
    /// Between a call to this operation and its adjoint, asserts that
    /// a given operation is called at most a certain number of times.
    ///
    /// Operation calls are considered, if they contain the the specified
    /// variant.  For example, if `op` is `X` also `Adjoint X` or `Controlled X`
    /// are counted, but if `op` is `Controlled X`, only `Controlled X`
    /// or `Controlled Adjoint X` are counted.
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
    /// ```qsharp
    /// within {
    ///     AllowAtMostNCallsCA(3, H, "Too many calls to H.");
    /// } apply {
    ///     use register = Qubit[4];
    ///     // Fails since this calls H four times, rather than the
    ///     // allowed maximum of three.
    ///     ApplyToEach(H, register);
    /// }
    /// ```
    ///
    /// # Remarks
    /// This operation may be replaced by a no-op on targets which do not
    /// support it.
    operation AllowAtMostNCallsCA<'TInput, 'TOutput>(
        nTimes : Int, op : ('TInput => 'TOutput),
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
    /// ```qsharp
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
