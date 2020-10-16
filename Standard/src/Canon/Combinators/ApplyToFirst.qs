// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    // # Overview
    // Variants of the operation that applies given one, two and three qubit
    // operation to the first one, two and three qubits of a register
    
    /// # Summary
    /// Applies an operation to the first qubit in the register.
    /// # Input
    /// ## op
    /// An operation to be applied to the first qubit
    /// ## register
    /// Qubit array to the first qubit of which the operation is applied
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubitA
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubitC
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubitCA
    operation ApplyToFirstQubit (op : (Qubit => Unit), register : Qubit[]) : Unit
    {
        if (Length(register) == 0)
        {
            fail $"Must have at least one qubit to act on.";
        }
        
        op(register[0]);
    }
    
    
    /// # Summary
    /// Applies an operation to the first qubit in the register.
    /// The modifier `A` indicates that the operation is adjointable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first qubit
    /// ## register
    /// Qubit array to the first qubit of which the operation is applied
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubit
    operation ApplyToFirstQubitA (op : (Qubit => Unit is Adj), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) == 0)
            {
                fail $"Must have at least one qubit to act on.";
            }
            
            op(register[0]);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Applies operation op to the first qubit in the register.
    /// The modifier `C` indicates that the operation is controllable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first qubit
    /// ## register
    /// Qubit array to the first qubit of which the operation is applied
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubit
    operation ApplyToFirstQubitC (op : (Qubit => Unit is Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) == 0)
            {
                fail $"Must have at least one qubit to act on.";
            }
            
            op(register[0]);
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Applies operation op to the first qubit in the register.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first qubit
    /// ## register
    /// Qubit array to the first qubit of which the operation is applied
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstQubit
    operation ApplyToFirstQubitCA (op : (Qubit => Unit is Adj + Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) == 0)
            {
                fail $"Must have at least one qubit to act on.";
            }
            
            op(register[0]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies an operation to the first two qubits in the register.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied to the first two qubits
    /// ## register
    /// Qubit array to the first two qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubitsA
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubitsC
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubitsCA
    operation ApplyToFirstTwoQubits (op : ((Qubit, Qubit) => Unit), register : Qubit[]) : Unit
    {
        if (Length(register) < 2)
        {
            fail $"Must have at least two qubits to act on.";
        }
        
        op(register[0], register[1]);
    }
    
    
    /// # Summary
    /// Applies an operation to the first two qubits in the register.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied to the first two qubits
    /// ## register
    /// Qubit array to the first two qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubits
    operation ApplyToFirstTwoQubitsA (op : ((Qubit, Qubit) => Unit is Adj), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 2)
            {
                fail $"Must have at least two qubits to act on.";
            }
            
            op(register[0], register[1]);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Applies an operation to the first two qubits in the register.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied to the first two qubits
    /// ## register
    /// Qubit array to the first two qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubits
    operation ApplyToFirstTwoQubitsC (op : ((Qubit, Qubit) => Unit is Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 2)
            {
                fail $"Must have at least two qubits to act on.";
            }
            
            op(register[0], register[1]);
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Applies an operation to the first two qubits in the register.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied to the first two qubits
    /// ## register
    /// Qubit array to the first two qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstTwoQubits
    operation ApplyToFirstTwoQubitsCA (op : ((Qubit, Qubit) => Unit is Adj + Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 2)
            {
                fail $"Must have at least two qubits to act on.";
            }
            
            op(register[0], register[1]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies an operation to the first three qubits in the register.
    /// # Input
    /// ## op
    /// An operation to be applied to the first three qubits
    /// ## register
    /// Qubit array to the first three qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1], register[2]);
    /// ```
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubitsC
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubitsA
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubitsCA
    operation ApplyToFirstThreeQubits (op : ((Qubit, Qubit, Qubit) => Unit), register : Qubit[]) : Unit
    {
        if (Length(register) < 3)
        {
            fail $"Must have at least three qubits to act on.";
        }
        
        op(register[0], register[1], register[2]);
    }
    
    
    /// # Summary
    /// Applies an operation to the first three qubits in the register.
    /// The modifier `A` indicates that the operation is adjointable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first three qubits
    /// ## register
    /// Qubit array to the first three qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1], register[2]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubits
    operation ApplyToFirstThreeQubitsA (op : ((Qubit, Qubit, Qubit) => Unit is Adj), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 3)
            {
                fail $"Must have at least three qubits to act on.";
            }
            
            op(register[0], register[1], register[2]);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Applies an operation to the first three qubits in the register.
    /// The modifier `C` indicates that the operation is controllable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first three qubits
    /// ## register
    /// Qubit array to the first three qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1], register[2]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubits
    operation ApplyToFirstThreeQubitsC (op : ((Qubit, Qubit, Qubit) => Unit is Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 3)
            {
                fail $"Must have at least three qubits to act on.";
            }
            
            op(register[0], register[1], register[2]);
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Applies an operation to the first three qubits in the register.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    /// # Input
    /// ## op
    /// An operation to be applied to the first three qubits
    /// ## register
    /// Qubit array to the first three qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1], register[2]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToFirstThreeQubits
    operation ApplyToFirstThreeQubitsCA (op : ((Qubit, Qubit, Qubit) => Unit is Adj + Ctl), register : Qubit[]) : Unit
    {
        body (...)
        {
            if (Length(register) < 3)
            {
                fail $"Must have at least three qubits to act on.";
            }
            
            op(register[0], register[1], register[2]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
}


