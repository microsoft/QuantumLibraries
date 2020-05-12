// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Given an array, returns a range over the indices of that array, suitable
    /// for use in a for loop.
    ///
    /// # Type Parameters
    /// ## 'TElement
    /// The type of elements of the array.
    ///
    /// # Input
    /// ## array
    /// An array for which a range of indices should be returned.
    ///
    /// # Output
    /// A range over all indices of the array.
    ///
    /// # Example
    /// The following `for` loops are equivalent:
    /// ```Q#
    /// for (idx in IndexRange(array)) { ... }
    /// for (idx in IndexRange(array)) { ... }
    /// ```
    function IndexRange<'TElement>(array : 'TElement[]) : Range {
       return 0..(Length(array) - 1);
    }

    internal function Identity<'T>(input : 'T) : 'T { return input; }

    /// # Summary
    /// Given an array, returns a new array containing elements of the original
    /// array along with the indices of each element.
    ///
    /// # Type Parameters
    /// ## 'TElement
    /// The type of elements of the array.
    ///
    /// # Input
    /// ## array
    /// An array whose elements are to be enumerated.
    ///
    /// # Output
    /// A new array containing elments of the original array along with their
    /// indices.
    ///
    /// # Example
    /// The following `for` loops are equivalent:
    /// ```Q#
    /// for (idx in IndexRange(array)) {
    ///     let element = array[idx];
    ///     ...
    /// }
    /// for ((idx, element) in Enumerated(array)) { ... }
    /// ```
    function Enumerated<'TElement>(array : 'TElement[]) : (Int, 'TElement)[] {
        return MappedByIndex(Identity<(Int, 'TElement)>, array);
    }

}
