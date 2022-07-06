// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Returns all consecutive subarrays of length `size`.
    ///
    /// # Description
    /// This function returns all `n - size + 1` subarrays of
    /// length `size` in order, where `n` is the length of `arr`.
    /// The first subarrays are `arr[0..size - 1], arr[1..size], arr[2..size + 1]`
    /// until the last subarray `arr[n - size..n - 1]`.
    ///
    /// If `size <= 0` or `size > n`, an empty array is returned.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    ///
    /// # Input
    /// ## size
    /// Length of the subarrays.
    ///
    /// ## array
    /// An array of elements.
    ///
    /// # Example
    /// ```qsharp
    /// // same as [[1, 2, 3], [2, 3, 4], [3, 4, 5]]
    /// let windows = Windows(3, [1, 2, 3, 4, 5]);
    /// ```
    function Windows<'T>(size : Int, array : 'T[]) : 'T[][] {
        let n = Length(array);

        if (size <= 0 or size > n) {
            return [];
        }

        mutable result = [[], size = n + 1 - size];

        for i in 0..n - size {
            set result w/= i <- array[i..i + size - 1];
        }

        return result;
    }

    /// # Summary
    /// Given an array, returns all its prefixes.
    ///
    /// # Description
    /// Returns an array of all prefixes, starting with an array that only
    /// has the first element until the complete array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    ///
    /// # Input
    /// ## array
    /// An array of elements.
    ///
    /// # Example
    /// ```qsharp
    /// let prefixes = Prefixes([23, 42, 144]);
    /// // prefixes = [[23], [23, 42], [23, 42, 144]]
    /// ```
    function Prefixes<'T>(array : 'T[]) : 'T[][] {
        return MappedOverRange(Prefix(_, array), IndexRange(array));
    }

    internal function Prefix<'T>(to : Int, array : 'T[]) : 'T[] {
        return array[0..to];
    }

    /// # Summary
    /// Applies an operation windowing over input but with constant target qubits.
    /// 
    /// # Input
    /// ## windowLen
    /// The size of each window
    /// ## op
    /// An operation with three arguments, one index, controls, and target. When applied it will be supplied with a windowed component of the `qubits` parameter, the `target` qubits (constant) and an index of the current window
    /// ## qubits
    /// The qubits the operation windows over
    /// ## target
    /// The target provided to each windowed operation
    operation ApplyWindowed(windowLen : Int, op : (Int, Qubit[], Qubit[]) => Unit,  qubits : Qubit[], target : Qubit[]) : Unit {
        let windows = Windows(windowLen, qubits); // Create windows of non-target qubits
        for (i, window) in Enumerated(windows) {
            op(i, window, target);
        }
    }

}
