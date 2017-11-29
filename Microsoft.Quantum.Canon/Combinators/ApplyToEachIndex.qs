// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Applies a single-qubit operation to each indexed element in a register.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    /// - @"microsoft.quantum.canon.ApplyToEachIndexA"
    /// - @"microsoft.quantum.canon.ApplyToEachIndexC"
    /// - @"microsoft.quantum.canon.ApplyToEachIndexCA"
    operation ApplyToEachIndex<'T>(singleElementOperation : ((Int, 'T) => ()), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..Length(register) - 1) {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.applytoeachindex"
    operation ApplyToEachIndexC<'T>(singleElementOperation : ((Int, 'T) => () : Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..Length(register) - 1) {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }

        controlled auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.applytoeachindex"
    operation ApplyToEachIndexA<'T>(singleElementOperation : ((Int, 'T) => () : Adjoint), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..Length(register) - 1) {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }

        adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.applytoeachindex"
    operation ApplyToEachIndexCA<'T>(singleElementOperation : ((Int, 'T) => () : Adjoint,Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..Length(register) - 1) {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

}
