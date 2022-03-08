// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies a single-qubit operation to each element in a register.
    /// The modifier `CA` indicates that the single-qubit operation is controllable
    /// and adjointable.
    ///
    /// # Input
    /// ## singleElementOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which the operation acts.
    ///
    /// # Example
    /// Prepare a three-qubit $\ket{+}$ state:
    /// ```qsharp
    /// using (register = Qubit[3]) {
    ///     ApplyToEachCA(H, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEach
    operation ApplyToEachCA<'T> (singleElementOperation : ('T => Unit is Adj + Ctl), register : 'T[])
    : Unit is Adj + Ctl {
        for idxQubit in IndexRange(register) {
            singleElementOperation(register[idxQubit]);
        }
    }


    /// # Summary
    /// Applies a single-qubit operation to each element in a register.
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
    /// The target on which the operation acts.
    ///
    /// # Example
    /// Prepare a three-qubit $\ket{+}$ state:
    /// ```qsharp
    /// using (register = Qubit[3]) {
    ///     ApplyToEachA(H, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEach
    operation ApplyToEachA<'T> (singleElementOperation : ('T => Unit is Adj), register : 'T[])
    : Unit is Adj {
        for idxQubit in IndexRange(register) {
            singleElementOperation(register[idxQubit]);
        }
    }


    /// # Summary
    /// Applies a single-qubit operation to each element in a register.
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
    /// The target on which the operation acts.
    ///
    /// # Example
    /// Prepare a three-qubit $\ket{+}$ state:
    /// ```qsharp
    /// using (register = Qubit[3]) {
    ///     ApplyToEachC(H, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEach
    operation ApplyToEachC<'T> (singleElementOperation : ('T => Unit is Ctl), register : 'T[])
    : Unit is Ctl {
        for idxQubit in IndexRange(register) {
            singleElementOperation(register[idxQubit]);
        }
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
    /// # Type Parameters
    /// ## 'T
    /// The target on which the operation acts.
    ///
    /// # Example
    /// Prepare a three-qubit $\ket{+}$ state:
    /// ```qsharp
    /// using (register = Qubit[3]) {
    ///     ApplyToEach(H, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToEachC
    /// - Microsoft.Quantum.Canon.ApplyToEachA
    /// - Microsoft.Quantum.Canon.ApplyToEachCA
    operation ApplyToEach<'T> (singleElementOperation : ('T => Unit), register : 'T[]) : Unit {
        for idxQubit in IndexRange(register) {
            singleElementOperation(register[idxQubit]);
        }
    }

}


