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
    ///
    /// let numbers = SequenceI(0, _); // function to create sequence from 0 to `to`
    /// let naturals = SequenceI(1, _); // function to create sequence from 1 to `to`
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
    /// The difference between `from` and `to` must fit into an `Int` value.
    ///
    /// ## Example
    /// ```qsharp
    /// let arr1 = SequenceL(0L, 3L); // [0L, 1L, 2L, 3L]
    /// let arr2 = SequenceL(23L, 29L); // [23L, 24L, 25L, 26L, 27L, 28L, 29L]
    /// let arr3 = SequenceL(-5L, -2L); // [-5L, -4L, -3L, -2L]
    /// ```
    function SequenceL (from : BigInt, to : BigInt) : BigInt[] {
        Fact(to >= from, $"`to` must be larger than `from`");
        Fact(to - from <= 0x07FFFFFFFFFFFFFFEL, $"different between `to` and `from` is too large");

        let length = BoolArrayAsInt(BigIntAsBoolArray(to - from)) + 1;
        mutable array = new BigInt[length];
        
        for (i in 0 .. length - 1) {
            set array w/= i <- from + IntAsBigInt(i);
        }

        return array;
    }

    function _RangeLength (range : Range) : Int {
        let start = RangeStart(range);
        let end = RangeEnd(range);
        let step = RangeStep(range);

        return ((end - start) / step) + 1;
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
    /// let arr4 = ArrayFromRange(-2..-2..-9); // [-2, -4, -6, -8]
    /// let arr5 = ArrayFromRange(-2..5..17); // [-2, 3, 8, 13]
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
