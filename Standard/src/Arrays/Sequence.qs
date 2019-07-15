// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Get an array of integers in a given interval.
    ///
    /// # Input
    /// ## from
    /// An inclusive nonnegative start index of the interval.
    /// ## to
    /// An inclusive nonnegative end index of the interval that is not smaller
    /// than `from`.
    ///
    /// # Output
    /// An array containing the sequence of numbers `from`, `from + 1`, ...,
    /// `to`.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let arr1 = SequenceI(0, 3); // [0, 1, 2, 3]
    /// let arr2 = SequenceI(23, 29); // [23, 24, 25, 26, 27, 28, 29]
    /// ```
    function SequenceI (from : Int, to : Int) : Int[] {
        let n = (to - from) + 1;
        mutable array = new Int[n];

        for (i in 0 .. n - 1) {
            set array w/= i <- from + i;
        }

        return array;
    }
}
