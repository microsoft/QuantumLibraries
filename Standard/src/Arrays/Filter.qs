// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Canon;

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
    ///
    /// # Example
    /// The following code demonstrates the "Filtered" function.
    /// A predicate is defined using the @"Microsoft.Quantum.Logical.GreaterThanI" function:
    /// ```qsharp
    /// open Microsoft.Quantum.Arrays;
    /// open Microsoft.Quantum.Logical;
    ///
    /// function FilteredDemo() : Unit {
    ///    let predicate = GreaterThanI(_, 5);
    ///    let filteredArray = Filtered(predicate, [2, 5, 9, 1, 8]);
    ///    Message($"{filteredArray}");
    /// }
    /// ```
    /// The outcome one should expect from this example will be an array of numbers greater than 5.
    function Filtered<'T>(predicate : ('T -> Bool), array : 'T[]) : 'T[] {
        mutable totalFound = 0;
        mutable idxArray = new Int[Length(array)];

        for idxElement in IndexRange(array) {
            if predicate(array[idxElement]) {
                set idxArray w/= totalFound <- idxElement;
                set totalFound += 1;
            }
        }

        return Subarray(idxArray[0 .. totalFound - 1], array);
    }

    /// # Summary
    /// Given a predicate and an array, returns the indices of that
    /// array where the predicate is true.
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
    /// An array of indices where `predicate` is true.
    function Where<'T>(predicate : ('T -> Bool), array : 'T[]) : Int[] {
        return Mapped(
            Fst,
            Filtered(
                Snd,
                Enumerated(Mapped(predicate, array))
            )
        );
    }

    /// # Summary
    /// Given an array and a predicate that is defined
    /// for the elements of the array, returns the number of elements
    /// an array that consists of those elements that satisfy the predicate.
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
    /// The number of elements in `array` that satisfy the predicate.
    ///
    /// # Example
    /// The following code demonstrates the "Count" function.
    /// A predicate is defined using the @"Microsoft.Quantum.Logical.GreaterThanI" function:
    /// ```qsharp
    ///  let predicate = GreaterThanI(_, 5);
    ///  let count = Count(predicate, [2, 5, 9, 1, 8]);
    ///  // count = 2
    /// ```
    function Count<'T>(predicate : ('T -> Bool), array : 'T[]) : Int {
        mutable totalFound = 0;

        for element in array {
            if predicate(element) {
                set totalFound += 1;
            }
        }

        return totalFound;
    }
}
