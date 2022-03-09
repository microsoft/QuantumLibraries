// Copyright (c) Microsoft Corporation.
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
    /// # Example
    /// Create three qubit state $\frac{1}{\sqrt{2}}\ket{0}\_2(\ket{0}\_1\ket{0}_3+\ket{1}\_1\ket{1}_3)$:
    /// ```qsharp
    ///     using (register = Qubit[3]) {
    ///         ApplyToSubregister(Exp([PauliX,PauliY],PI() / 4.0,_), [1,3], register);
    ///     }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToSubregisterA
    /// - Microsoft.Quantum.Canon.ApplyToSubregisterC
    /// - Microsoft.Quantum.Canon.ApplyToSubregisterCA
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
    /// - Microsoft.Quantum.Canon.ApplyToSubregister
    operation ApplyToSubregisterA (op : (Qubit[] => Unit is Adj), idxs : Int[], target : Qubit[]) : Unit
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
    /// - Microsoft.Quantum.Canon.ApplyToSubregister
    operation ApplyToSubregisterC (op : (Qubit[] => Unit is Ctl), idxs : Int[], target : Qubit[]) : Unit
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
    /// - Microsoft.Quantum.Canon.ApplyToSubregister
    operation ApplyToSubregisterCA (op : (Qubit[] => Unit is Ctl + Adj), idxs : Int[], target : Qubit[]) : Unit
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregister
    function RestrictedToSubregisterA (op : (Qubit[] => Unit is Adj), idxs : Int[]) : (Qubit[] => Unit is Adj)
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregister
    function RestrictedToSubregisterC (op : (Qubit[] => Unit is Ctl), idxs : Int[]) : (Qubit[] => Unit is Ctl)
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RestrictedToSubregister
    function RestrictedToSubregisterCA (op : (Qubit[] => Unit is Adj + Ctl), idxs : Int[]) : (Qubit[] => Unit is Adj + Ctl)
    {
        return ApplyToSubregisterCA(op, idxs, _);
    }
    
}


