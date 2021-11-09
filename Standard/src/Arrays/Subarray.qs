// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Takes an array and a list of locations and
    /// produces a new array formed from the elements of the original
    /// array that match the given locations.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have
    /// an array `'T[]` and a list of locations `Int[]` defining the subarray.
    /// The construction of the subarray is a based on generating a new, deep
    /// copy of the given array as opposed to maintaining references.
    ///
    /// If `indices` contains repeated elements, the corresponding elements 
    /// of `array` will likewise be repeated.
    /// If all elements of `indices` are unique, this function will return 
    /// a subset of `array` if `Length(indices) < Length(array)`, or
    /// a permutation of `array` if `indices` and `array` are of the same length.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    ///
    /// # Input
    /// ## indices
    /// A list of integers that is used to define the subarray.
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `out` of elements whose indices correspond to the subarray,
    /// such that `out[idx] == array[indices[idx]]`.
    /// 
    /// # Example
    /// 
    /// ```qsharp
    /// open Microsoft.Quantum.Arrays;
    ///
    /// function SubarrayDemo() : Unit {
    ///     let array = [1, 2, 3, 4];
    ///     let permutation = Subarray([3, 0, 2, 1], array); // [4, 1, 3, 2]
    ///     let duplicates = Subarray([1, 2, 2], array);     // [2, 3, 3]
    /// }
    /// ```
    function Subarray<'T> (indices : Int[], array : 'T[]) : 'T[] {
        let nSliced = Length(indices);
        mutable sliced = new 'T[nSliced];

        for idx in 0 .. nSliced - 1 {
            set sliced w/= idx <- array[indices[idx]];
        }

        return sliced;
    }

}
