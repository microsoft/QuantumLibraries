// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies a single-qubit operation to each indexed element in a register.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEach
    /// - Microsoft.Quantum.Canon.ApplyToEachIndexA
    /// - Microsoft.Quantum.Canon.ApplyToEachIndexC
    /// - Microsoft.Quantum.Canon.ApplyToEachIndexCA
    operation ApplyToEachIndex<'T> (singleElementOperation : (Int, 'T) => Unit, register : 'T[])
    : Unit {
        for idxQubit in IndexRange(register) {
            singleElementOperation(idxQubit, register[idxQubit]);
        }
    }


    /// # Summary
    /// Applies a single-qubit operation to each indexed element in a register.
    /// The modifier `C` indicates that the single-qubit operation is controllable.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEachIndex
    operation ApplyToEachIndexC<'T> (singleElementOperation : (Int, 'T) => Unit is Ctl, register : 'T[])
    : Unit is Ctl {
        for idxQubit in IndexRange(register) {
            singleElementOperation(idxQubit, register[idxQubit]);
        }
    }


    /// # Summary
    /// Applies a single-qubit operation to each indexed element in a register.
    /// The modifier `A` indicates that the single-qubit operation is adjointable.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEachIndex
    operation ApplyToEachIndexA<'T> (singleElementOperation : ((Int, 'T) => Unit is Adj), register : 'T[])
    : Unit is Adj {
        for idxQubit in IndexRange(register) {
            singleElementOperation(idxQubit, register[idxQubit]);
        }
    }


    /// # Summary
    /// Applies a single-qubit operation to each indexed element in a register.
    /// The modifier `CA` indicates that the single-qubit operation is adjointable
    /// and controllable.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEachIndex
    operation ApplyToEachIndexCA<'T> (singleElementOperation : ((Int, 'T) => Unit is Adj + Ctl), register : 'T[])
    : Unit is Adj + Ctl {
        for idxQubit in IndexRange(register) {
            singleElementOperation(idxQubit, register[idxQubit]);
        }
    }

}

