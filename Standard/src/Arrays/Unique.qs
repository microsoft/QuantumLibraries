// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    /// # Summary
    /// Returns a new array that has no equal adjacent elements.
    ///
    /// # Description
    /// Given some array of elements and a function to test equality, this
    /// function returns a new array in which the relative order of elements
    /// is kept, but all adjacent elements which are equal are filtered to
    /// just a single element.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## equal
    /// A function that compares two elements such that `a` is considered to
    /// be equal to `b` if `equal(a, b)` is `true`.
    /// ## array
    /// The array to be filtered for unique elements.
    ///
    /// # Output
    /// Array with no equal adjacent elements.
    ///
    /// # Remarks
    /// If there are multiple elements that are equal but not next to each other,
    /// there will be multiple occurrences in the output array.  Use this function
    /// together with `Sorted` to get an array with overall unique elements.
    ///
    /// # Example
    /// ```qsharp
    /// let unique1 = Unique(EqualI, [1, 1, 3, 3, 2, 5, 5, 5, 7]);
    /// // same as [1, 3, 2, 5, 7]
    /// let unique2 = Unique(EqualI, [2, 2, 1, 1, 2, 2, 1, 1]);
    /// // same as [2, 1, 2, 1];
    /// let unique3 = Unique(EqualI, Sorted(LessThanOrEqualI, [2, 2, 1, 1, 2, 2, 1, 1]));
    /// // same as [1, 2];
    /// ```
    function Unique<'T>(equal : (('T, 'T) -> Bool), array : 'T[]) : 'T[] {
        if Length(array) == 0 {
            return [];
        }

        mutable unique = ConstantArray(Length(array), Head(array));
        mutable count = 1;

        for elem in Rest(array) {
            if (not equal(elem, unique[count - 1])) {
                set unique w/= count <- elem;
                set count += 1;
            }
        }

        return unique[0..count - 1];
    }
}
