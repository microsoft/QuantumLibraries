// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Given an array, returns a Range over the indices of that array, suitable
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
    /// for (idx in IndexRange(arr)) { ... }
    /// for (idx in IndexRange(arr)) { ... }
    /// ```
    function IndexRange<'TElement>(array : 'TElement[]) : Range {
       return 0..(Length(array) - 1);
    }

}
