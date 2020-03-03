// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Given two arrays of the same type and a predicate that is defined
    /// for pairs of elements of the arrays, checks whether the arrays are equal.
    ///
    /// # Remarks
    /// This function is defined for generic types; i.e., whenever we have
    /// two arrays `'T[]` and a function `equal: ('T, 'T) -> Bool`, this function returns
    /// a `Bool` value that indicates whether the arrays are equal.
    /// For two arrays to be considered equal, they have to have equal length
    /// and the elements in the same positions in the arrays have to be equal.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each array's elements.
    ///
    /// # Input
    /// ## predicate
    /// A function from tuple `('T, 'T)` to `Bool` that is used to check whether two elements of arrays are equal.
    /// ## array1
    /// The first array to be compared.
    /// ## array2
    /// The second array to be compared.
    ///
    /// # Output
    /// The value `true` if and only if `array1` and `array2` are equal.
    /// That is, if both arrays have the same length, and all elements are equal
    /// as defined by `equal`.
    ///
    /// # Example 
    /// The following code checks whether different pairs of arrays are equal:
    /// ```qsharp
    /// open Microsoft.Quantum.Arrays;
    /// open Microsoft.Quantum.Logical;
    ///
    /// function EqualADemo() : Unit {
    ///     let equalArrays = EqualA(EqualI, [2, 3, 4], [2, 3, 4]);
    ///     let differentLength = EqualA(EqualD, [2.0, 3.0, 4.0], [2.0, 3.0]);
    ///     let differentElements = EqualA(EqualR, [One, Zero], [One, One]);
    ///     Message($"Equal arrays are {equalArrays ? "equal" | "not equal"}");
    ///     Message($"Arrays of different length are {differentLength ? "equal" | "not equal"}");
    ///     Message($"Arrays of the same length with different elements are {differentElements ? "equal" | "not equal"}");
    /// }
    /// ```
    ///
    function EqualA<'T>(equal : (('T, 'T) -> Bool), array1 : 'T[], array2 : 'T[]) : Bool {
        if (Length(array1) != Length(array2)) {
            return false;
        }
        return All(equal, Zip(array1, array2));
    }

}
