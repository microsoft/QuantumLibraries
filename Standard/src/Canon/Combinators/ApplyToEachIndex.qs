// Copyright (c) Microsoft Corporation. All rights reserved.
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
    operation ApplyToEachIndex<'T> (singleElementOperation : ((Int, 'T) => Unit), register : 'T[]) : Unit
    {
        for (idxQubit in IndexRange(register))
        {
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
    operation ApplyToEachIndexC<'T> (singleElementOperation : ((Int, 'T) => Unit : Controlled), register : 'T[]) : Unit
    {
        body (...)
        {
            for (idxQubit in IndexRange(register))
            {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }
        
        controlled distribute;
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
    operation ApplyToEachIndexA<'T> (singleElementOperation : ((Int, 'T) => Unit : Adjoint), register : 'T[]) : Unit
    {
        body (...)
        {
            for (idxQubit in IndexRange(register))
            {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }
        
        adjoint invert;
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
    operation ApplyToEachIndexCA<'T> (singleElementOperation : ((Int, 'T) => Unit : Adjoint, Controlled), register : 'T[]) : Unit
    {
        body (...)
        {
            for (idxQubit in IndexRange(register))
            {
                singleElementOperation(idxQubit, register[idxQubit]);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
}


