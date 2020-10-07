// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Given two arrays, returns a new array of pairs such that each pair
    /// contains an element from each original array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the left array elements.
    /// ## 'U
    /// The type of the right array elements.
    ///
    /// # Input
    /// ## left
    /// An array containing values for the first element of each tuple.
    /// ## right
    /// An array containing values for the second element of each tuple.
    ///
    /// # Output
    /// An array containing pairs of the form `(left[idx], right[idx])` for
    /// each `idx`. If the two arrays are not of equal length, the output will
    /// be as long as the shorter of the inputs.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let left = [1, 3, 71];
    /// let right = [false, true];
    /// let pairs = Zipped(left, right); // [(1, false), (3, true)]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Zipped3
    /// - Microsoft.Quantum.Arrays.Zipped4
    /// - Microsoft.Quantum.Arrays.Unzipped
    function Zipped<'T, 'U>(left : 'T[], right : 'U[]) : ('T, 'U)[] {
        let nElements = Length(left) < Length(right)
                        ? Length(left)
                        | Length(right);
        mutable output = new ('T, 'U)[nElements];

        for (idxElement in 0 .. nElements - 1) {
            set output w/= idxElement <- (left[idxElement], right[idxElement]);
        }

        return output;
    }

    /// # Summary
    /// Given three arrays, returns a new array of 3-tuples such that each 3-tuple
    /// contains an element from each original array.
    ///
    /// # Type Parameters
    /// ## 'T1
    /// The type of the first array elements.
    /// ## 'T2
    /// The type of the second array elements.
    /// ## 'T3
    /// The type of the third array elements.
    ///
    /// # Input
    /// ## first
    /// An array containing values for the first element of each tuple.
    /// ## second
    /// An array containing values for the second element of each tuple.
    /// ## third
    /// An array containing values for the third element of each tuple.
    ///
    /// # Output
    /// An array containing 3-tuples of the form `(first[idx], second[idx], third[idx])` for
    /// each `idx`. If the three arrays are not of equal length, the output will
    /// be as long as the shorter of the inputs.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Zipped
    /// - Microsoft.Quantum.Arrays.Zipped4
    function Zipped3<'T1, 'T2, 'T3> (first : 'T1[], second : 'T2[], third : 'T3[]) : ('T1, 'T2, 'T3)[] {
        let nElements = Min([Length(first), Length(second), Length(third)]);
        mutable output = new ('T1, 'T2, 'T3)[nElements];

        for (idxElement in 0 .. nElements - 1) {
            set output w/= idxElement <- (first[idxElement], second[idxElement], third[idxElement]);
        }

        return output;
    }

    /// # Summary
    /// Given four arrays, returns a new array of 4-tuples such that each 4-tuple
    /// contains an element from each original array.
    ///
    /// # Type Parameters
    /// ## 'T1
    /// The type of the first array elements.
    /// ## 'T2
    /// The type of the second array elements.
    /// ## 'T3
    /// The type of the third array elements.
    /// ## 'T4
    /// The type of the fourth array elements.
    ///
    /// # Input
    /// ## first
    /// An array containing values for the first element of each tuple.
    /// ## second
    /// An array containing values for the second element of each tuple.
    /// ## third
    /// An array containing values for the third element of each tuple.
    /// ## fourth
    /// An array containing values for the fourth element of each tuple.
    ///
    /// # Output
    /// An array containing 4-tuples of the form `(first[idx], second[idx], third[idx], fourth[idx])` for
    /// each `idx`. If the four arrays are not of equal length, the output will
    /// be as long as the shorter of the inputs.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Zipped
    /// - Microsoft.Quantum.Arrays.Zipped3
    function Zipped4<'T1, 'T2, 'T3, 'T4> (first : 'T1[], second : 'T2[], third : 'T3[], fourth : 'T4[]) : ('T1, 'T2, 'T3, 'T4)[] {
        let nElements = Min([Length(first), Length(second), Length(third), Length(fourth)]);
        mutable output = new ('T1, 'T2, 'T3, 'T4)[nElements];

        for (idxElement in 0 .. nElements - 1) {
            set output w/= idxElement <- (first[idxElement], second[idxElement], third[idxElement], fourth[idxElement]);
        }

        return output;
    }

    /// # Summary
    /// Given an array of 2-tuples, returns a tuple of two arrays, each containing
    /// the elements of the tuples of the input array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first element in each tuple
    /// ## 'U
    /// The type of the second element in each tuple
    ///
    /// # Input
    /// ## arr
    /// An array containing 2-tuples
    ///
    /// # Output
    /// Two arrays, the first one containing all first elements of the input
    /// tuples, the second one containing all second elements of the input tuples.
    ///
    /// # Example
    /// ```Q#
    /// // split is same as ([6, 5, 5, 3, 2, 1], [true, false, false, false, true, false])
    /// let split = Unzipped([(6, true), (5, false), (5, false), (3, false), (2, true), (1, false)]);
    /// ```
    ///
    /// # Remark
    /// This function is equivalent to `(Mapped(Fst<'T, 'U>, arr), Mapped(Snd<'T, 'U>, arr))`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Zipped
    function Unzipped<'T, 'U>(arr : ('T, 'U)[]) : ('T[], 'U[]) {
        let nElements = Length(arr);
        mutable first = new 'T[nElements];
        mutable second = new 'U[nElements];
        for (idxElement in 0 .. nElements - 1) {
            let (left, right) = arr[idxElement];
            set first w/= idxElement <- left;
            set second w/= idxElement <- right;
        }
        return (first, second);
    }
}


