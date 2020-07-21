// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Internal function used to fail with meaningful error messages.
    ///
    /// # Remarks
    /// This function is intended to be emulated by simulation runtimes, so as
    /// to allow forwarding formatted contradictions.
    internal function FormattedFailure<'T>(actual : 'T, expected : 'T, message : String) : Unit {
        fail $"{message}\n\tExpected:\t{expected}\n\tActual:\t{actual}";
    }

}
