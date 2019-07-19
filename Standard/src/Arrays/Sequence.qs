// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Get an array of integers in a given interval.
    ///
    /// # Input
    /// ## from
    /// An inclusive start index of the interval.
    /// ## to
    /// An inclusive end index of the interval that is not smaller than `from`.
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
    /// let arr3 = SequenceI(-5, -2); // [-5, -4, -3, -2]
    /// ```
    function SequenceI (from : Int, to : Int) : Int[] {
        Fact(to >= from, $"`to` must be larger than `from`");

        let n = (to - from) + 1;
        mutable array = new Int[n];

        for (i in 0 .. n - 1) {
            set array w/= i <- from + i;
        }

        return array;
    }

    /// # Summary
    /// Get an array of integers in an interval specified with start and length.
    ///
    /// # Input
    /// ## from
    /// An inclusive nonnegative start index of the interval.
    /// ## length
    /// A nonnegative integer specifying the length of the resulting array.
    ///
    /// # Output
    /// An array containing the sequence of numbers `from`, `from + 1`, ...,
    /// `from + length - 1`.  If `length` is 0, the resulting array is empty.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let arr1 = SequenceL(0L, 4); // [0L, 1L, 2L, 3L]
    /// let arr2 = SequenceL(23L, 7); // [23L, 24L, 25L, 26L, 27L, 28L, 29L]
    /// let arr3 = SequenceL(-5L, 4); // [-5L, -4L, -3L, -2L]
    /// let arr4 = SequenceL(10L, 0); // []
    /// ```
    function SequenceL (from : BigInt, length : Int) : BigInt[] {
        Fact(length >= 0, $"`length` must be nonnegative");

        mutable array = new BigInt[length];
        
        for (i in 0 .. length - 1) {
            set array w/= i <- from + IntAsBigInt(i);
        }

        return array;
    }

    /// # Summary
    /// Get a sequence of numbers starting from 0.
    ///
    /// # Input
    /// ## count
    /// A nonnegative number of how many elements the resulting array should
    /// contain.
    ///
    /// # Output
    /// An array with the elements `0`, `1`, ..., `count - 1`.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// let arr1 = Numbers(3); // [0, 1, 2]
    /// let arr2 = Numbers(5); // [0, 1, 2, 3, 4]
    /// let arr3 = Numbers(0); // []
    /// ```
    function Numbers (count : Int) : Int[] {
        Fact(count >= 0, $"`count` must be nonnegative");

        if (count == 0) {
            return new Int[0];
        } else {
            return SequenceI(0, count - 1);
        }
    }

    function _RangeLength (range : Range) : Int {
        mutable length = 0;
        for (i in range) {
            set length = length + 1;
        }
        return length;
    }

    /// # Summary
    /// Given a range, returns an array containing the numbers visited in the
    /// range.
    ///
    /// # Input
    /// ## range
    /// A Q# range
    ///
    /// # Output
    /// An array containing each number visited in the range in the same order.
    ///
    /// # Remarks
    /// This function needs to iterate through the range twice.  For simple
    /// ranges with an increment of `1`, it's faster to use `SequenceI`.
    ///
    /// ## Example
    /// ```qsharp
    /// let arr1 = ArrayFromRange(1..5); // [1, 2, 3, 4, 5]
    /// let arr2 = ArrayFromRange(5..-1..1); // [5, 4, 3, 2, 1]
    /// let arr3 = ArrayFromRange(13..2..19); // [13, 15, 17, 19]
    /// ```
    function ArrayFromRange (range : Range) : Int[] {
        mutable array = new Int[_RangeLength(range)];
        mutable i = 0;

        for (elem in range) {
            set array w/= i <- elem;
            set i = i + 1;
        }

        return array;
    }
}
