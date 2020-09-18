// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Returns all consecutive subarrays of length `size`.
    ///
    /// # Description
    /// This function returns all `n - size + 1` subarrays of
    /// length `size` in order, where `n` is the length of `arr`.
    /// The first subarrays are `arr[0..size - 1], arr[1..size], arr[2..size + 1]`
    /// until the last subarray `arr[n - size..n - 1]`.
    ///
    /// If `size <= 0` or `size > n`, an empty array is returned.
    ///
    /// # Input
    /// ## size
    /// Length of the subarrays
    ///
    /// ## arr
    /// Array
    ///
    /// # Example
    /// ```Q#
    /// // same as [[1, 2, 3], [2, 3, 4], [3, 4, 5]]
    /// let windows = Windows(3, [1, 2, 3, 4, 5]);
    /// ```
    function Windows<'T>(size : Int, arr : 'T[]) : 'T[][] {
        let n = Length(arr);

        if (size <= 0 or size > n) {
            return new 'T[][0];
        }

        mutable result = new 'T[][n + 1 - size];

        for (i in 0..n - size) {
            set result w/= i <- arr[i..i + size - 1];
        }

        return result;
    }
}
