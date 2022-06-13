// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Logical;

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
        return array[1 ...];
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
        return array[... Length(array) - 2];
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
        return ElementAt(_, array);
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
    /// Returns a tuple of first and all remaining elements of the array.
    ///
    /// # Type Parameters
    /// ## 'A
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// An array with at least one element.
    ///
    /// # Output
    /// A tuple of first and all remaining elements of the array.
    function HeadAndRest<'A>(array : 'A[]) : ('A, 'A[]) {
        return (Head(array), Rest(array));
    }

    /// # Summary
    /// Returns a tuple of all but one and the last element of the array.
    ///
    /// # Type Parameters
    /// ## 'A
    /// The type of the array elements.
    ///
    /// # Input
    /// ## array
    /// An array with at least one element.
    ///
    /// # Output
    /// A tuple of all but one and the last element of the array.
    function MostAndTail<'A>(array : 'A[]) : ('A[], 'A) {
        return (Most(array), Tail(array));
    }

    /// # Summary
    /// Creates an array of given length with all elements equal to given value.
    ///
    /// # Input
    /// ## length
    /// Length of the new array.
    /// ## value
    /// The value of each element of the new array.
    ///
    /// # Output
    /// A new array of length `length`, such that every element is `value`.
    ///
    /// # Example
    /// The following code creates an array of 3 Boolean values, each equal to `true`:
    /// ```qsharp
    /// let array = ConstantArray(3, true);
    /// ```
    function ConstantArray<'T> (length : Int, value : 'T) : 'T[] {
        return [value, size = length];
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
    /// # Example
    /// ```qsharp
    /// let array = [10, 11, 12, 13, 14, 15];
    /// // The following line returns [10, 12, 15].
    /// let subarray = Excluding([1, 3, 4], array);
    /// ```
    function Excluding<'T>(remove : Int[], array : 'T[]) : 'T[] {
        let nSliced = Length(remove);
        let nElements = Length(array);

        if nElements - nSliced <= 0 {
            return [];
        }

        //Would be better with sort function
        //Or way to add elements to array
        mutable arrayKeep = SequenceI(0, nElements - 1);
        mutable sliced = [array[0], size = nElements - nSliced];
        mutable counter = 0;

        for idx in remove {
            set arrayKeep w/= idx <- -1;
        }

        for idx in arrayKeep {
            if idx >= 0 {
                set sliced w/= counter <- array[idx];
                set counter += 1;
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
    /// # Example
    /// ```qsharp
    /// let array = [10, 11, 12];
    /// // The following line returns [10, 12, 15, 2, 2, 2].
    /// let output = Padded(-6, 2, array);
    /// // The following line returns [2, 2, 2, 10, 12, 15].
    /// let output = Padded(6, 2, array);
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
    /// Splits an array into multiple parts of equal length.
    ///
    /// # Input
    /// ## nElements
    /// The length of each chunk. Must be positive.
    /// ## arr
    /// The array to be split.
    ///
    /// # Output
    /// A array containing each chunk of the original array.
    ///
    /// # Remarks
    /// Note that the last element of the output may be shorter
    /// than `nElements` if `Length(arr)` is not divisible by `nElements`.
    function Chunks<'T>(nElements : Int, arr : 'T[]) : 'T[][] {
        Fact(nElements > 0, "nElements must be positive");
        mutable output = [];
        mutable remaining = arr;
        while (not IsEmpty(remaining)) {
            let nElementsToTake = MinI(Length(remaining), nElements);
            set output += [remaining[...nElementsToTake - 1]];
            set remaining = remaining[nElementsToTake...];
        }
        return output;
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
    /// will contain all remaining elements. If the array is split exactly, the
    /// last array will be the empty array, indicating there are no remaining elements.
    /// In other words, `Tail(Partitioned(...))` will always return the remaining
    /// elements, while `Most(Partitioned(...))` will always return the complete
    /// partitions of the array.
    ///
    /// # Example
    /// ```qsharp
    /// // The following returns [[1, 5], [3], [7]];
    /// let split = Partitioned([2,1], [1,5,3,7]);
    /// // The following returns [[1, 5], [3, 7], []];
    /// let split = Partitioned([2,2], [1,5,3,7]);
    /// ```
    function Partitioned<'T>(nElements: Int[], arr: 'T[]) : 'T[][] {
        mutable output = [[], size = Length(nElements) + 1];
        mutable currIdx = 0;
        for idx in IndexRange(nElements) {
            if currIdx + nElements[idx] > Length(arr) {
                fail "Partitioned argument out of bounds.";
            }
            set output w/= idx <- arr[currIdx..currIdx + nElements[idx] - 1];
            set currIdx = currIdx + nElements[idx];
        }
        set output w/= Length(nElements) <- arr[currIdx..Length(arr) - 1];
        return output;
    }

    /// # Summary
    /// Returns true if and only if an array is empty.
    ///
    /// # Input
    /// ## array
    /// The array to be checked.
    ///
    /// # Output
    /// `true` if and only if the array is empty (has length 0).
    function IsEmpty<'T>(array : 'T[]) : Bool {
        return Length(array) == 0;
    }

    internal function IsValuePresent(permutation : Int[], value : Int) : Bool {
        let index = IndexOf(EqualI(value, _), permutation);
        return index != -1;
    }

    /// # Summary
    /// Outputs true if and only if a given array represents a permutation.
    ///
    /// # Description
    /// Given an array `array` of length `n`, returns true if and only if
    /// each integer from `0` to `n - 1` appears exactly once in `array`, such
    /// that `array` can be interpreted as a permutation on `n` elements.
    ///
    /// # Input
    /// ## permutation
    /// An array that may or may not represent a permutation.
    ///
    /// # Output
    /// `true` if and only if the array is a permutation.
    ///
    /// # Remarks
    /// An array of length zero is trivially a permutation.
    ///
    /// # Example
    /// The following Q# code prints the message "All diagnostics completed
    /// successfully":
    /// ```qsharp
    /// Fact(IsPermutation([2, 0, 1], "");
    /// Contradiction(IsPermutation([5, 0, 1], "[5, 0, 1] isn't a permutation");
    /// Message("All diagnostics completed successfully.");
    /// ```
    function IsPermutation(permutation : Int[]) : Bool {
        return All(IsValuePresent(permutation, _), RangeAsIntArray(IndexRange(permutation)));
    }

    // NB: This function is internal, but not marked as internal so as to allow
    //     unit tests to check its behavior. In the future, tests should be
    //     redesigned to check only publicly accessible behavior.
    /// # Summary
    /// Returns the order elements in an array need to be swapped to produce an ordered array.
    /// Assumes swaps occur in place.
    ///
    /// # Input
    /// ## newOrder
    /// Array with the permutation of the indices of the new array. There should be $n$ elements,
    /// each being a unique integer from $0$ to $n-1$.
    ///
    /// # Output
    /// The tuple represents the two indices to be swapped. The swaps begin at the lowest index.
    ///
    /// # Example
    /// ```qsharp
    /// // The following returns [(0, 5),(0, 4),(0, 1),(0, 3)];
    /// let swapOrder = _SwapOrderToPermuteArray([5, 3, 2, 0, 1, 4]);
    /// ```
    ///
    /// # Remarks
    /// ## Pseudocode
    /// ```
    /// for index in 0..Length(newOrder) - 1 {
    ///     while newOrder[index] != index {
    ///         Switch newOrder[index] with newOrder[newOrder[index]]
    ///     }
    /// }
    /// ```
    function _SwapOrderToPermuteArray(newOrder : Int[]) : (Int, Int)[] {
        // Check to verify the new ordering actually is a permutation of the indices
        Fact(IsPermutation(newOrder), $"The new ordering is not a permutation");

        mutable swaps = [];
        mutable order = newOrder;

        // for each value, whenever the index and value don't match, swap until it does
        for index in IndexRange(order) {
            while (not EqualI(order[index], index)) {
                set swaps += [(index, order[index])];
                set order = Swapped(order[index], index, order);
            }
        }

        return swaps;
    }

    /// # Summary
    /// Applies a swap of two elements in an array.
    ///
    /// # Input
    /// ## firstIndex
    /// Index of the first element to be swapped.
    ///
    /// ## secondIndex
    /// Index of the second element to be swapped.
    ///
    /// ## arr
    /// Array with elements to be swapped.
    ///
    /// # Output
    /// The array with the in place swap applied.
    ///
    /// # Example
    /// ```qsharp
    /// // The following returns [0, 3, 2, 1, 4]
    /// Swapped(1, 3, [0, 1, 2, 3, 4]);
    /// ```
    function Swapped<'T>(firstIndex : Int, secondIndex : Int, arr : 'T[]) : 'T[] {
        Fact(firstIndex >= 0 and firstIndex < Length(arr), "First index is out of bounds");
        Fact(secondIndex >= 0 and secondIndex < Length(arr), "Second index is out of bounds");
        return arr
            w/ firstIndex <- arr[secondIndex]
            w/ secondIndex <- arr[firstIndex];
    }

    /// # Summary
    /// Turns a list of 2-tuples into a nested array.
    ///
    /// # Input
    /// ## tupleList
    /// List of 2-tuples to be turned into a nested array.
    ///
    /// # Output
    /// A nested array with length matching the tupleList.
    ///
    /// # Example
    /// ```qsharp
    /// // The following returns [[2, 3], [4, 5]]
    /// TupleArrayAsNestedArray([(2, 3), (4, 5)]);
    /// ```
    function TupleArrayAsNestedArray<'T>(tupleList : ('T, 'T)[]) : 'T[][] {
        mutable newArray = [[], size = Length(tupleList)];
        for idx in IndexRange(tupleList) {
            let (tupleLeft, tupleRight) = tupleList[idx];
            set newArray w/= idx <- [tupleLeft, tupleRight];
        }
        return newArray;
    }

}
