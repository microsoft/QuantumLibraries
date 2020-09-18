// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Interleaves two arrays of (almost) same size.
    ///
    /// # Description
    /// This function returns the interleaving of two arrays, starting
    /// with the first element from the first array, then the first
    /// element from the second array, and so on.
    ///
    /// The first array must either be
    /// of the same length as the second one, or can have one more element.
    ///
    /// # Input
    /// ## first
    /// The first array to be interleaved.
    ///
    /// ## second
    /// The second array to be interleaved.
    ///
    /// # Output
    /// Interleaved array
    ///
    /// # Example
    /// ```Q#
    /// // same as int1 = [1, -1, 2, -2, 3, -3]
    /// let int1 = Interleaved([1, 2, 3], [-1, -2, -3])
    ///
    /// // same as int2 = [false, true, false, true, false]
    /// let int2 = Interleaved(ConstantArray(3, false), ConstantArray(2, true));
    /// ```
    function Interleaved<'T>(first : 'T[], second : 'T[]) : 'T[] {
        let lFirst = Length(first);
        let lSecond = Length(second);

        Fact(lFirst >= lSecond and lFirst - lSecond <= 1, "Array `first` is either of same size as `second`, or has one more element");

        return new 'T[lFirst + lSecond]
            w/ 0..2..(lFirst + lSecond - 1) <- first
            w/ 1..2..(lFirst + lSecond - 1) <- second;
    }
}
