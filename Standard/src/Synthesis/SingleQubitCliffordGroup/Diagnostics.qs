namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Requires that a single-qubit Clifford operator is equal to the identity.
    ///
    /// # Input
    /// ## actual
    /// The value to be checked.
    /// ## message
    /// Failure message string to be used when `actual` is not equal to the identity.
    function IdentityFact1C(op : SingleQubitClifford, message : String) : Unit {
        EqualityFact1C(op, Identity1C(), message);
    }

    /// # Summary
    /// Requires that a single-qubit Clifford operator has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The value to be checked.
    /// ## expected
    /// The expected value.
    /// ## message
    /// Failure message string to be used when `actual` is not equal to `expected`.
    function EqualityFact1C(actual : SingleQubitClifford, expected : SingleQubitClifford, message : String)
    : Unit {
        if actual::E != expected::E or actual::S != expected::S or actual::X != expected::X or actual::Omega != expected::Omega {
            FormattedFailure(actual, expected, message);
        }
    }

}
