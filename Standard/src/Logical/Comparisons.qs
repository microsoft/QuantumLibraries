// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Logical {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Given a comparison function, returns a new function that
    /// lexographically compares two arrays.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the elements of the arrays being compared.
    ///
    /// # Input
    /// ## elementComparison
    /// A function that compares two elements `x` and `y` and returns if
    /// `x` is less than or equal to `y`.
    ///
    /// # Output
    /// A function that compares two arrays `xs` and `ys` and returns if
    /// `xs` occurs before or equal to `ys` in lexographical ordering.
    ///
    /// # Remarks
    /// The lexographic comparison between two arrays `xs` and `ys` is defined
    /// by the following procedure. We say that two elements `x` and `y`
    /// are equivalent if `elementComparison(x, y)` and `elementComparison(y, x)`
    /// are both true.
    ///
    /// - Both arrays are compared element-by-element until the first pair of
    ///   elements that are not equivalent. The array containing the element
    ///   that occurs first according to `elementComparison` is said to occur
    ///   first in lexographical ordering.
    /// - If no inequivalent elements are found, and one array is longer than
    ///   the other, the shorter array is said to occur first.
    ///
    /// # Examples
    /// ```Q#
    /// let arrayComparison = LexographicComparison(LessThanOrEqualD);
    /// let data = [
    ///     [1.1, 2.2, 3.3],
    ///     [1.1, 2.2],
    ///     [0.2, 2.2],
    ///     [1.1, 2.7]
    /// ];
    /// let sorted = Sorted(arrayComparison, data);
    /// // sorted:
    /// // [
    /// //     [0.2, 2.2],
    /// //     [1.1, 2.2],
    /// //     [1.1, 2.2, 3.3],
    /// //     [1.1, 2.7]
    /// // ];
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Sorted
    function LexographicComparison<'T>(elementComparison : (('T, 'T) -> Bool)) : (('T[], 'T[]) -> Bool) {
        return LessThanLexographic(elementComparison, _, _);
    }

    /// # Summary
    /// Used to implement `LexographicComparison`.
    internal function LessThanLexographic<'T>(comparison : (('T, 'T) -> Bool), left : 'T[], right : 'T[]) : Bool {
        for ((l, r) in Zip(left, right)) {
            let lessThanOrEqual = comparison(l, r);
            let greaterThanOrEqual = comparison(r, l);
            let equal = lessThanOrEqual and greaterThanOrEqual;
            if (lessThanOrEqual and not equal) {
                return true;
            } elif (greaterThanOrEqual and not equal) {
                return false;
            }
        }

        // At this point, all items in the common prefix of both arrays
        // are equal to each other under comparison (l ≤ r and r ≤ l).
        // Thus, if left is shorter than or equal to right, then left occurs
        // at or before right in lexographical ordering.
        return Length(left) <= Length(right);
    }

}
