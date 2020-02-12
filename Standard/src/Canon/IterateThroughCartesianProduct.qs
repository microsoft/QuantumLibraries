// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies an operation for each index in the Cartesian product of several
    /// ranges.
    ///
    /// # Description
    /// Iteratively applies an operation for each element of the Cartesian product
    /// of `0..(bounds[0] - 1)`, `0..(bounds[1] - 1)`, ..., `0..(bounds[Length(bounds) - 1] - 1)`
    ///
    /// # Input
    /// ## bounds
    /// An array specifying the ranges to be iterated over, with each range
    /// being specified as an integer length.
    /// ## op
    /// An operation to be called for each element of the given Cartesian product.
    ///
    /// # Example
    /// Given an operation `op`, the following two snippets are equivalent:
    /// ```Q#
    /// IterateThroughCartesianProduct([3, 4, 5], op);
    /// ```
    /// ```Q#
    /// op([0, 0, 0]);
    /// op([1, 0, 0]);
    /// op([2, 0, 0]);
    /// op([0, 1, 0]);
    /// // ...
    /// op([0, 3, 0]);
    /// op([0, 0, 1]);
    /// //
    /// op([2, 3, 4]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.IterateThroughCartesianPower
    operation IterateThroughCartesianProduct(bounds : Int[], op : (Int[] => Unit)) : Unit {
        mutable arr = new Int[Length(bounds)];
        mutable finished = false;

        repeat {
            if (not finished) {
                op(arr);
            }
        }
        until (finished)
        fixup {
            //computes the next element in the Cartesian product
            set arr w/= 0 <- arr[0] + 1;

            for (i in 0 .. Length(arr) - 2) {
                if (arr[i] == bounds[i]) {
                    set arr w/= i + 1 <- arr[i + 1] + 1;
                    set arr w/= i <- 0;
                }
            }

            if (arr[Length(arr) - 1] == bounds[Length(arr) - 1]) {
                set finished = true;
            }
        }
    }


    /// # Summary
    /// Applies an operation for each index in the Cartesian power of an
    /// integer range.
    ///
    /// # Description
    /// Iteratively applies an operation for each element of a Cartesian power
    /// of the range `0..(bound - 1)`.
    ///
    /// # Input
    /// ## power
    /// The Cartesian power to which the range `0..(bound - 1)` should be
    /// raised.
    /// ## bound
    /// A specification of the range to be iterated over, given as the length
    /// of the range.
    /// ## op
    /// An operation to be called for each element of the given Cartesian power.
    ///
    /// # Example
    /// Given an operation `op`, the following two snippets are equivalent:
    /// ```Q#
    /// IterateThroughCartesianPower(2, 3, op);
    /// ```
    /// ```Q#
    /// op([0, 0]);
    /// op([1, 0]);
    /// op([2, 0]);
    /// op([0, 1]);
    /// // ..
    /// op([2, 2]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.IterateThroughCartesianProduct
    operation IterateThroughCartesianPower (power : Int, bound : Int, op : (Int[] => Unit)) : Unit {
        IterateThroughCartesianProduct(ConstantArray(power, bound), op);
    }

}
