// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies an operation to a subregister of a register, with qubits
    /// specified by an array of their indices.
    ///
    /// # Input
    /// ## op
    /// Operation to apply to subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # Remarks
    /// ## Example
    /// Create three qubit state $\frac{1}{\sqrt{2}}\ket{0}\_2(\ket{0}\_1\ket{0}_3+\ket{1}\_1\ket{1}_3)$:
    /// ```qsharp
    ///     using (register = Qubit[3]) {
    ///         ApplyToSubregister(Exp([PauliX,PauliY],PI() / 4.0,_), [1,3], register);
    ///     }
    /// ```
    ///
    /// # See Also
    /// - ApplyToSubregisterA
    /// - ApplyToSubregisterC
    /// - ApplyToSubregisterCA
    operation ApplyToSubregister (op : (Qubit[] => Unit), idxs : Int[], target : Qubit[]) : Unit {
        let subregister = Subarray(idxs, target);
        op(subregister);
    }


    /// # Summary
    /// Applies an operation to a subregister of a register, with qubits
    /// specified by an array of their indices.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to apply to subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterA (op : (Qubit[] => Unit : Adjoint), idxs : Int[], target : Qubit[]) : Unit
    {
        body (...)
        {
            ApplyToSubregister(op, idxs, target);
        }
        
        adjoint (...)
        {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
    }
    
    
    /// # Summary
    /// Applies an operation to a subregister of a register, with qubits
    /// specified by an array of their indices.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// Operation to apply to subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterC (op : (Qubit[] => Unit : Controlled), idxs : Int[], target : Qubit[]) : Unit
    {
        body (...)
        {
            ApplyToSubregister(op, idxs, target);
        }
        
        controlled (controls, ...)
        {
            let cop = Controlled op;
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }
    
    
    /// # Summary
    /// Applies an operation to a subregister of a register, with qubits
    /// specified by an array of their indices.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to apply to subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterCA (op : (Qubit[] => Unit : Controlled, Adjoint), idxs : Int[], target : Qubit[]) : Unit
    {
        body (...)
        {
            ApplyToSubregister(op, idxs, target);
        }
        
        adjoint (...)
        {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
        
        controlled (controls, ...)
        {
            let cop = Controlled op;
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
        
        controlled adjoint (controls, ...)
        {
            let cop = Controlled (Adjoint op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }
    
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.
    ///
    /// # Input
    /// ## op
    /// Operation to be restricted to a subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be restricted.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregisterA
    /// - Microsoft.Quantum.Canon.RestrictedToSubregisterC
    /// - Microsoft.Quantum.Canon.RestrictedToSubregisterCA
    function RestrictedToSubregister (op : (Qubit[] => Unit), idxs : Int[]) : (Qubit[] => Unit) {
        return ApplyToSubregister(op, idxs, _);
    }
    
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to be restricted to a subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be restricted.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - RestrictedToSubregister
    function RestrictedToSubregisterA (op : (Qubit[] => Unit : Adjoint), idxs : Int[]) : (Qubit[] => Unit : Adjoint)
    {
        return ApplyToSubregisterA(op, idxs, _);
    }
    
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// Operation to be restricted to a subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be restricted.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregister
    function RestrictedToSubregisterC (op : (Qubit[] => Unit : Controlled), idxs : Int[]) : (Qubit[] => Unit : Controlled)
    {
        return ApplyToSubregisterC(op, idxs, _);
    }
    
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to be restricted to a subregister.
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be restricted.
    /// ## target
    /// Register on which the operation acts.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregister
    function RestrictedToSubregisterCA (op : (Qubit[] => Unit : Adjoint, Controlled), idxs : Int[]) : (Qubit[] => Unit : Adjoint, Controlled)
    {
        return ApplyToSubregisterCA(op, idxs, _);
    }
    
}


