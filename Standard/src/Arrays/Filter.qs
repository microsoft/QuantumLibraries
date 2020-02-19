// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Given an array and a predicate that is defined
    /// for the elements of the array, returns an array that consists of
    /// those elements that satisfy the predicate.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have
    /// an array `'T[]` and a predicate `'T -> Bool` we can filter elements.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    ///
    /// # Input
    /// ## predicate
    /// A function from `'T` to Boolean that is used to filter elements.
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `'T[]` of elements that satisfy the predicate.
    function Filtered<'T> (predicate : ('T -> Bool), array : 'T[]) : 'T[] {
        mutable totalFound = 0;
        mutable idxArray = new Int[Length(array)];

        for (idxElement in IndexRange(array)) {
            if (predicate(array[idxElement])) {
                set idxArray w/= totalFound <- idxElement;
                set totalFound = totalFound + 1;
            }
        }

        return Subarray(idxArray[0 .. totalFound - 1], array);
    }


    /// # Example
    /// The following code demonstrates the "Filtered" function.
    /// A predicate is defined using the @"microsoft.quantum.logical.greaterthani" function:
    /// open Microsoft.Quantum.Arrays;
    /// open Microsoft.Quantum.Logical;
    ///
    /// function FilteredDemo () : Unit {
    ///    let predicate = GreaterThanI(_, 5);
    ///    let filteredArray = Filtered(predicate, [2, 5, 9, 1, 8]);
    ///    Message($"{filteredArray}");
    /// }

    /// ```

}
