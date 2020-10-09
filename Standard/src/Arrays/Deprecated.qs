// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Logical;

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
    /// let pairs = Zip(left, right); // [(1, false), (3, true)]
    /// ```
    ///
    /// # See Also
    /// - Zip3
    /// - Zip4
    /// - Unzipped
    @Deprecated("Microsoft.Quantum.Arrays.Zipped")
    function Zip<'T, 'U> (left : 'T[], right : 'U[]) : ('T, 'U)[] {
        return Zipped(left, right);
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
    /// - Zip
    /// - Zip4
    @Deprecated("Microsoft.Quantum.Arrays.Zipped3")
    function Zip3<'T1, 'T2, 'T3> (first : 'T1[], second : 'T2[], third : 'T3[]) : ('T1, 'T2, 'T3)[] {
        return Zipped3(first, second, third);
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
    /// - Zip
    /// - Zip3
    @Deprecated("Microsoft.Quantum.Arrays.Zipped4")
    function Zip4<'T1, 'T2, 'T3, 'T4> (first : 'T1[], second : 'T2[], third : 'T3[], fourth : 'T4[]) : ('T1, 'T2, 'T3, 'T4)[] {
        return Zipped4(first, second, third, fourth);
    }

    

    /// # Summary
    /// Returns an array containing the elements of another array,
    /// excluding elements at a given list of indices.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the array elements.
    ///
    /// # Input
    /// ## remove
    /// An array of indices denoting which elements should be excluded
    /// from the output.
    /// ## array
    /// Array of which the values in the output array are taken.
    ///
    /// # Output
    /// An array `output` such that `output[0]` is the first element
    /// of `array` whose index does not appear in `remove`,
    /// such that `output[1]` is the second such element, and so
    /// forth.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let array = [10, 11, 12, 13, 14, 15];
    /// // The following line returns [10, 12, 15].
    /// let subarray = Exclude([1, 3, 4], array);
    /// ```
    function Exclude<'T>(remove : Int[], array : 'T[]) : 'T[] {
        return Excluding(remove, array);
    }

}
