// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Returns the representative of an integer on the ring
    /// Z / nZ, where the representative is guaranteed to be in the
    /// interval [0, n).
    ///
    /// # Remarks
    /// This differs from the `%` operator in how negative `x` are handled.
    internal function RingRepresentative(x : Int, n : Int) : Int {
        Fact(n > 0, "Modulus must be positive.");
        return (x + n) % n;
    }

}
