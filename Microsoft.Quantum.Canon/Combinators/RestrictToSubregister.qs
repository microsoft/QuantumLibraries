// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    
    /// # Summary
    /// Applies an operation to an array of indices of a register, i.e., a subregister.  
    ///
    /// # Input
    /// ## op
    /// Operation to be applied to a subregister. 
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied. 
    /// ## target 
    /// Register on which the operation acts. 
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregistera"
    /// - @"microsoft.quantum.canon.applytosubregisterc"
    /// - @"microsoft.quantum.canon.applytosubregisterca"
    operation ApplyToSubregister(op : (Qubit[] => ()), idxs : Int[], target : Qubit[]) : () {
        body {
            let subregister = Subarray(idxs, target);
            op(subregister);
        }
    }

    /// # Summary
    /// Applies an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'A' indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to be applied to a subregister. 
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied. 
    /// ## target 
    /// Register on which the operation acts. 
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
    }

    /// # Summary
    /// Applies an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'C' indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// Operation to be applied to a subregister. 
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied. 
    /// ## target 
    /// Register on which the operation acts. 
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }

    /// # Summary
    /// Applies an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'CA' indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// Operation to be applied to a subregister. 
    /// ## idxs
    /// Array of indices, indicating to which qubits the operation will be applied. 
    /// ## target 
    /// Register on which the operation acts. 
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytosubregister"
    operation ApplyToSubregisterCA(op : (Qubit[] => () : Controlled, Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
        controlled adjoint (controls) {
            let cop = (Controlled Adjoint op);
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
    /// - @"microsoft.quantum.canon.restricttosubregistera"
    /// - @"microsoft.quantum.canon.restricttosubregisterc"
    /// - @"microsoft.quantum.canon.restricttosubregisterca"
    function RestrictToSubregister(op : (Qubit[] => ()), idxs : Int[]) : (Qubit[] => ()) {
        return ApplyToSubregister(op, idxs, _);
    }
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'A' indicates that the operation is adjointable.
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
    /// - @"microsoft.quantum.canon.restricttosubregister"
    function RestrictToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[]) : (Qubit[] => () : Adjoint) {
        return ApplyToSubregisterA(op, idxs, _);
    }
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'C' indicates that the operation is controllable.
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
    /// - @"microsoft.quantum.canon.restricttosubregister"
    function RestrictToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[]) : (Qubit[] => () : Controlled) {
        return ApplyToSubregisterC(op, idxs, _);
    }
    
    /// # Summary
    /// Restricts an operation to an array of indices of a register, i.e., a subregister.  
    /// The modifier 'CA' indicates that the operation is controllable and adjointable.
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
    /// - @"microsoft.quantum.canon.restricttosubregister"
    function RestrictToSubregisterCA(op : (Qubit[] => () : Adjoint, Controlled), idxs : Int[]) : (Qubit[] => () : Adjoint, Controlled) {
        return ApplyToSubregisterCA(op, idxs, _);
    }

}
