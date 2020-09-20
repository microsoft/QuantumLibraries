// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Given an array and a function that is defined
    /// for the elements of the array, returns a new array that consists
    /// of the images of the original array under the function.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have
    /// an array `'T[]` and a function `mapper: 'T -> 'U` we can map the elements
    /// of the array and produce a new array of type `'U[]`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    /// ## 'U
    /// The result type of the `mapper` function.
    ///
    /// # Input
    /// ## mapper
    /// A function from `'T` to `'U` that is used to map elements.
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `'U[]` of elements that are mapped by the `mapper` function.
    function Mapped<'T, 'U> (mapper : ('T -> 'U), array : 'T[]) : 'U[] {
        mutable resultArray = new 'U[Length(array)];

        for (idxElement in IndexRange(array)) {
            set resultArray w/= idxElement <- mapper(array[idxElement]);
        }

        return resultArray;
    }

    /// # Summary
    /// Given an array and a function that is defined
    /// for the indexed elements of the array, returns a new array that consists
    /// of the images of the original array under the function.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have
    /// an array `'T[]` and a function `mapper: (Int, 'T) -> 'U` we can map the elements
    /// of the array and produce a new array of type `'U[]`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    /// ## 'U
    /// The result type of the `mapper` function.
    ///
    /// # Input
    /// ## mapper
    /// A function from `(Int, 'T)` to `'U` that is used to map elements
    /// and their indices.
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `'U[]` of elements that are mapped by the `mapper` function.
    ///
    /// # Remarks
    /// ## Example
    /// The following two lines are equivalent:
    /// ```qsharp
    /// let arr = MapIndex(f, [x0, x1, x2]);
    /// ```
	/// and
	/// ```qsharp
    /// let arr = [f(0, x0), f(1, x1), f(2, x2)];
    /// ```
    ///
    /// # See Also
    /// - Mapped
    function MappedByIndex<'T, 'U> (mapper : ((Int, 'T) -> 'U), array : 'T[]) : 'U[] {
        mutable resultArray = new 'U[Length(array)];

        for (idxElement in IndexRange(array)) {
            set resultArray w/= idxElement <- mapper(idxElement, array[idxElement]);
        }

        return resultArray;
    }

    /// # Summary
    /// Given a range and a function that takes an integer as input,
    /// returns a new array that consists
    /// of the images of the range values under the function.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have
    /// a function `mapper: Int -> 'T` we can map the values
    /// of the range and produce an array of type `'T[]`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The result type of the `mapper` function.
    ///
    /// # Input
    /// ## mapper
    /// A function from `Int` to `'T` that is used to map range values.
    /// ## range
    /// A range of integers.
    ///
    /// # Output
    /// An array `'T[]` of elements that are mapped by the `mapper` function.
    function MappedOverRange<'T> (mapper : (Int -> 'T), range : Range) : 'T[] {
        let start = RangeStart(range);
        let step = RangeStep(range);
        let end = RangeEnd(range);
        if ((end - start) / step >= 0) {
            let nTerms = (end - start) / step + 1;
            mutable resultArray = new 'T[nTerms];
            mutable idxElement = 0;
            for (elem in range) {
                set resultArray w/= idxElement <- mapper(elem);
                set idxElement += 1;
            }
            return resultArray;
        } else {
            return new 'T[0];
        }
    }

    /// # Summary
    /// Given an array and an operation that is defined
    /// for the elements of the array, returns a new array that consists
    /// of the images of the original array under the operation.
    ///
    /// # Remarks
    /// The operation is defined for generic types, i.e., whenever we have
    /// an array `'T[]` and an operation `action : 'T -> 'U` we can map the elements
    /// of the array and produce a new array of type `'U[]`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    /// ## 'U
    /// The result type of the `action` operation.
    ///
    /// # Input
    /// ## action
    /// An operation from `'T` to `'U` that is applied to each element.
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `'U[]` of elements that are mapped by the `action` operation.
    operation ForEach<'T, 'U> (action : ('T => 'U), array : 'T[]) : 'U[] {
        mutable resultArray = new 'U[Length(array)];

        for (idxElement in IndexRange(array)) {
            set resultArray w/= idxElement <- action(array[idxElement]);
        }

        return resultArray;
    }

}

