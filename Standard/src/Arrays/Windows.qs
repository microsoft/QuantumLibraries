// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Canon;

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
    /// Applies an operation windowing over an input register.
    /// 
    /// # Input
    /// ## windowLen
    /// The size of each window.
    /// ## op
    /// An operation on a register that will be provided with the current window and its index.
    /// ## register
    /// The register the operation windows over.
    ///
    /// # Example
    /// The example below shows how to use `ApplyToEachWindow` to construct a parity function
    /// ```qsharp
    /// operation Parity(qubits : Qubit[], target : Qubit) : Unit {
    ///     ApplyToEachWindow(1, (_, q) => CNOT(q[0], target), qubits);
    /// }
    /// ```
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of register elements.
    operation ApplyToEachWindow<'T>(windowLen : Int, op : (Int, 'T[]) => Unit, register : 'T[]) : Unit {
        ApplyToEach(op, Enumerated(Windows(windowLen, register)));
    }

    /// # Summary
    /// Applies an operation windowing over an input register. The modifier `A` indicates that the single-qubit operation is adjointable.
    /// 
    /// # Input
    /// ## windowLen
    /// The size of each window.
    /// ## op
    /// An operation on a register that will be provided with the current window and its index.
    /// ## register
    /// The register the operation windows over.
    ////
    /// # Type Parameters
    /// ## 'T
    /// The type of register elements.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.ApplyToEachWindow
    operation ApplyToEachWindowA<'T>(windowLen : Int, op : (Int, 'T[]) => Unit is Adj, register : 'T[]) : Unit is Adj {
        ApplyToEachA(op, Enumerated(Windows(windowLen, register)));
    }

    /// # Summary
    /// Applies an operation windowing over an input register. The modifier `C` indicates that the single-qubit operation is controllable.
    /// 
    /// # Input
    /// ## windowLen
    /// The size of each window.
    /// ## op
    /// An operation on a register that will be provided with the current window and its index.
    /// ## register
    /// The register the operation windows over.
    ////
    /// # Type Parameters
    /// ## 'T
    /// The type of register elements.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.ApplyToEachWindow
    operation ApplyToEachWindowC<'T>(windowLen : Int, op : (Int, 'T[]) => Unit is Ctl, register : 'T[]) : Unit is Ctl {
        ApplyToEachC(op, Enumerated(Windows(windowLen, register)));
    }

    /// # Summary
    /// Applies an operation windowing over an input register. The modifier `CA` indicates that the single-qubit operation is controllable and adjointable.
    /// 
    /// # Input
    /// ## windowLen
    /// The size of each window.
    /// ## op
    /// An operation on a register that will be provided with the current window and its index.
    /// ## register
    /// The register the operation windows over.
    ////
    /// # Type Parameters
    /// ## 'T
    /// The type of register elements.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.ApplyToEachWindow
    operation ApplyToEachWindowCA<'T>(windowLen : Int, op : (Int, 'T[]) => Unit is Adj + Ctl, register : 'T[]) : Unit is Adj + Ctl {
        ApplyToEachCA(op, Enumerated(Windows(windowLen, register)));
    }
}
