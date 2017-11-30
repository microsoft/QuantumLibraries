// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachCA<'T>(singleElementOperation : ('T => () : Adjoint, Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleElementOperation(register[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachA<'T>(singleElementOperation : ('T => ():Adjoint), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleElementOperation(register[idxQubit]);
            }
        }

        adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachC<'T>(singleElementOperation : ('T => ():Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleElementOperation(register[idxQubit]);
            }
        }

        controlled auto
    }

    /// # Summary
    /// Applies a single-qubit operation to each element in a register.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Remarks
    /// ## Example
    /// Prepare a three-qubit $\ket{+} state:
    /// ```qsharp
    ///     using (register = Qubit[3]) {
    ///         ApplyToEach(H, register);
    ///     }
    /// ```
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEachA"
    /// - @"microsoft.quantum.canon.ApplyToEachC"
    /// - @"microsoft.quantum.canon.ApplyToEachCA"
    operation ApplyToEach<'T>(singleElementOperation : ('T => ()), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleElementOperation(register[idxQubit]);
            }
        }
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    function ApplyToEachF<'T>( func : ('T -> () ), array : 'T[])  : ()
    {
        for ( idx in 0..(Length(array) - 1)) {
            func(array[idx]);
        }
    }

}
