// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Given two sorted arrays, returns a single array containing the
    /// elements of both in sorted order. Used internally by merge sort.
    internal function Merged<'T>(comparison : (('T, 'T) -> Bool), left : 'T[], right : 'T[]) : 'T[] {
        mutable result = new 'T[0];
        mutable l = left;
        mutable r = right;
        while ((not IsEmpty(l)) and (not IsEmpty(r))) {
            if (comparison(Head(l), Head(r))) {
                set result += [Head(l)];
                set l = Rest(l);
            } else {
                set result += [Head(r)];
                set r = Rest(r);
            }
        }

        // Note that at this point, either or both of l and r are empty,
        // such that we can simply append both to our result to get the
        // whole merged array.
        return result + l + r;
    }

    /// # Summary
    /// Given an array, returns whether that array is sorted as defined by
    /// a given comparison function.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## comparison
    /// A function that compares two elements such that `a` is considered to
    /// be less than or equal to `b` if `comparison(a, b)` is `true`.
    /// ## array
    /// The array to be checked.
    ///
    /// # Output
    /// `true` if and only if for each pair of elements `a` and `b` of
    /// `array` occurring in that order, `comparison(a, b)` is `true`.
    ///
    /// # Remarks
    /// The function `comparison` is assumed to be transitive, such that
    /// if `comparison(a, b)` and `comparison(b, c)`, then `comparison(a, c)`
    /// is assumed. If this property does not hold, then the output of this
    /// function may be incorrect.
    function IsSorted<'T>(comparison : (('T, 'T) -> Bool), array : 'T[]) : Bool {
        return All(
            comparison,
            Zipped(Most(array), Rest(array))
        );
    }

    /// # Summary
    /// Given an array, returns the elements of that array sorted by a given
    /// comparison function.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## comparison
    /// A function that compares two elements such that `a` is considered to
    /// be less than or equal to `b` if `comparison(a, b)` is `true`.
    /// ## array
    /// The array to be sorted.
    ///
    /// # Output
    /// An array containing the same elements as `array`, such that for all
    /// elements `a` occurring earlier than elements `b`, `comparison(a, b)`
    /// is `true`.
    ///
    /// # Example
    /// The following snippet sorts an array of integers to occur in ascending
    /// order:
    /// ```qsharp
    /// let sortedArray = Sorted(LessThanOrEqualI, [3, 17, 11, -201, -11]);
    /// ```
    /// 
    /// # Remarks
    /// The function `comparison` is assumed to be transitive, such that
    /// if `comparison(a, b)` and `comparison(b, c)`, then `comparison(a, c)`
    /// is assumed. If this property does not hold, then the output of this
    /// function may be incorrect.
    ///
    /// As this is a function, the results are completely deterministic, even
    /// when two elements are considered equal under `comparison`;
    /// that is, when `comparison(a, b)` and `comparison(b, a)` are both `true`.
    /// In particular, the sort performed by this function is guaranteed to be
    /// stable, so that if two elements `a` and `b` occur in that order within
    /// `array` and are considered equal under `comparison`, then `a` will also
    /// appear before `b` in the output.
    ///
    /// For example:
    /// ```qsharp
    /// function LastDigitLessThanOrEqual(left : Int, right : Int) : Bool {
    ///     return LessThanOrEqualI(
    ///         left % 10, right % 10
    ///     );
    /// }
    ///
    /// function SortedByLastDigit() : Int[] {
    ///     return Sorted(LastDigitLessThanOrEqual, [3, 37, 11, 17]);
    /// }
    /// // returns [11, 3, 37, 17].
    /// ```
    function Sorted<'T>(comparison : (('T, 'T) -> Bool), array : 'T[]) : 'T[] {
        if Length(array) <= 1 {
            return array;
        } else {
            let idxPivot = Length(array) / 2;
            let left = array[...idxPivot - 1];
            let right = array[idxPivot...];

            // Sort each sublist, then merge them back into a single combined
            // list and return.
            return Merged(
                comparison, 
                Sorted(comparison, left),
                Sorted(comparison, right)
            );
        }
    }

}
