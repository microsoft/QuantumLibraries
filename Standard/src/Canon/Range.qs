// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Returns true if and only if input range is empty.
    ///
    /// # Input
    /// ## rng
    /// Any range
    ///
    /// # Output
    /// True, if and only if `rng` is empty
    ///
    /// # Remark
    /// This function needs to check at most one range index
    /// to determine whether the range is empty.
    function IsRangeEmpty(rng : Range) : Bool {
        for (idx in rng) {
            return false;
        }
        return true;
    }
}
