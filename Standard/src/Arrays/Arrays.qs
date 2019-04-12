// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Create an array that contains the same elements as an input array but in Reversed
    /// order.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// An array whose elements are to be copied in Reversed order.
    ///
    /// # Output
    /// An array containing the elements `array[Length(array) - 1]` .. `array[0]`.
    function Reversed<'T>(array : 'T[]) : 'T[] {
        let nElements = Length(array);
        return array[nElements-1..-1..0];
    }

    /// # Summary
    /// Creates an array that is equal to an input array except that the first array
    /// element is dropped.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// An array whose second to last elements are to form the output array.
    ///
    /// # Output
    /// An array containing the elements `array[1..Length(array) - 1]`.
    function Rest<'T> (array : 'T[]) : 'T[] {
        return array[1 .. Length(array) - 1];
    }

    /// # Summary
    /// Creates an array that is equal to an input array except that the last array
    /// element is dropped.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// An array whose first to second-to-last elements are to form the output array.
    ///
    /// # Output
    /// An array containing the elements `array[0..Length(array) - 2]`.
    function Most<'T> (array : 'T[]) : 'T[] {
        return array[0 .. Length(array) - 2];
    }

    function LookupImpl<'T> (array : 'T[], index : Int) : 'T {
        return array[index];
    }

    /// # Summary
    /// Given an array, returns a function which returns elements of that
    /// array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the elements of the array being represented as a lookup
    /// function.
    ///
    /// # Input
    /// ## array
    /// The array to be represented as a lookup function.
    ///
    /// # Output
    /// A function `f` such that `f(idx) == f[idx]`.
    ///
    /// # Remarks
    /// This function is mainly useful for interoperating with functions and
    /// operations that take `Int -> 'T` functions as their arguments. This
    /// is common, for instance, in the generator representation library,
    /// where functions are used to avoid the need to record an entire array
    /// in memory.
    function LookupFunction<'T> (array : 'T[]) : (Int -> 'T) {
        return LookupImpl(array, _);
    }

    /// # Summary
    /// Returns the last element of the array.
    ///
    /// # Type Parameters
    /// ## 'A
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// Array of which the last element is taken. Array must have length at least 1.
    ///
    /// # Output
    /// The last element of the array.
    function Tail<'A> (array : 'A[]) : 'A {
        EqualityFactB(Length(array) > 0, true, $"Array must be of the length at least 1");
        return array[Length(array) - 1];
    }

    /// # Summary
    /// Returns the first element of the array.
    ///
    /// # Type Parameters
    /// ## 'A
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// Array of which the first element is taken. Array must have length at least 1.
    ///
    /// # Output
    /// The first element of the array.
    function Head<'A> (array : 'A[]) : 'A {
        EqualityFactB(Length(array) > 0, true, $"Array must be of the length at least 1");
        return array[0];
    }

    /// # Summary
    /// Creates an array of given length with all elements equal to given value.
    ///
    /// # Input
    /// ## length
    /// Length of the new array.
    /// ## value
    /// A value that will be contained at each index of the new array.
    ///
    /// # Output
    /// A new array of length `length`, such that every element is `value`.
    function ConstantArray<'T> (length : Int, value : 'T) : 'T[] {
        mutable arr = new 'T[length];

        for (i in 0 .. length - 1) {
            set arr[i] = value;
        }

        return arr;
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
    function Exclude<'T> (remove : Int[], array : 'T[]) : 'T[] {
        let nSliced = Length(remove);
        let nElements = Length(array);

        //Would be better with sort function
        //Or way to add elements to array
        mutable arrayKeep = new Int[nElements];
        mutable sliced = new 'T[nElements - nSliced];
        mutable counter = 0;

        for (idx in 0 .. nElements - 1) {
            set arrayKeep[idx] = idx;
        }

        for (idx in 0 .. nSliced - 1) {
            set arrayKeep[remove[idx]] = -1;
        }

        for (idx in 0 .. nElements - 1) {
            if (arrayKeep[idx] >= 0) {
                set sliced[counter] = array[arrayKeep[idx]];
                set counter = counter + 1;
            }
        }

        return sliced;
    }

    /// # Summary
    /// Returns an array padded at with specified values up to a
    /// specified length.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the array elements.
    ///
    /// # Input
    /// ## nElementsTotal
    /// The length of the padded array. If this is positive, `inputArray`
    /// is padded at the head. If this is negative, `inputArray` is padded
    /// at the tail.
    /// ## defaultElement
    /// Default value to use for padding elements.
    /// ## inputArray
    /// Array whose values are at the head of the output array.
    ///
    /// # Output
    /// An array `output` that is the `inputArray` padded at the head
    /// with `defaultElement`s until `output` has length `nElementsTotal`
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let array = [10, 11, 12];
    /// // The following line returns [10, 12, 15, 2, 2, 2].
    /// let output = Padded(-6, array, 2);
    /// // The following line returns [2, 2, 2, 10, 12, 15].
    /// let output = Padded(6, array, 2);
    /// ```
    function Padded<'T> (nElementsTotal : Int, defaultElement : 'T, inputArray : 'T[]) : 'T[] {
        let nElementsInitial = Length(inputArray);
        let nAbsElementsTotal = AbsI(nElementsTotal);
        EqualityFactB(nAbsElementsTotal >= nElementsInitial, true, $"Specified output array length must be longer than `inputArray` length.");
        let nElementsPad = nAbsElementsTotal - nElementsInitial;
        let padArray = ConstantArray(nElementsPad, defaultElement);

        return nElementsTotal >= 0
               ? padArray + inputArray  // Padded at head.
               | inputArray + padArray; // Padded at tail.
    }

    /// # Summary
    /// Splits an array into multiple parts.
    ///
    /// # Input
    /// ## nElements
    /// Number of elements in each part of array
    /// ## arr
    /// Input array to be split.
    ///
    /// # Output
    /// Multiple arrays where the first array is the first `nElements[0]` of `arr`
    /// and the second array are the next `nElements[1]` of `arr` etc. The last array
    /// will contain all remaining elements.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// // The following returns [[1,5],[3],[7]];
    /// let split = SplitArray([2,1], [1,5,3,7]);
    /// ```
    function Partitioned<'T>(nElements: Int[], arr: 'T[]) : 'T[][] {
        mutable output = new 'T[][Length(nElements)+1];
        mutable currIdx = 0;
        for (idx in IndexRange(nElements)) {
            if(currIdx + nElements[idx] > Length(arr)) {
                fail "SplitArray argument out of bounds.";
            }
            set output[idx] = arr[currIdx..currIdx + nElements[idx]-1];
            set currIdx = currIdx + nElements[idx];
        }
        set output[Length(nElements)] = arr[currIdx..Length(arr)-1];
        return output;
    }

}


