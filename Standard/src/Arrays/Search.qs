// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Returns the first index of the first element in an array that satisfies
    /// a given predicate. If no such element exists, returns -1.
    ///
    /// # Input
    /// ## predicate
    /// A predicate function acting on elements of the array.
    /// ## arr
    /// An array to be searched using the given predicate.
    ///
    /// # Output
    /// Either the smallest index `idx` such that `predicate(arr[idx])` is true,
    /// or -1 if no such element exists.
    ///
    /// # Example
    /// Suppose that `IsEven : Int -> Bool` is a function that returns `true`
    /// if and only if its input is even. Then, this can be used with `IndexOf`
    /// to find the first even element in an array:
    /// ```qsharp
    /// let items = [1, 3, 17, 2, 21];
    /// let idx = IndexOf(IsEven, items); // returns 3
    /// ```
    function IndexOf<'T>(predicate : ('T -> Bool), arr : 'T[]) : Int {
        for (idx in IndexRange(arr)) {
            if (predicate(arr[idx])) {
                return idx;
            }
        }
        return -1;
    }

}
