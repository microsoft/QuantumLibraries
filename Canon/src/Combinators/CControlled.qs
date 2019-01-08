// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    /// # Summary
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is true. If false, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfC
    /// - Microsoft.Quantum.Canon.ApplyIfA
    /// - Microsoft.Quantum.Canon.ApplyIfCA
    operation ApplyIf<'T> (op : ('T => Unit), bit : Bool, target : 'T) : Unit
    {
        if (bit)
        {
            op(target);
        }
    }
    
    
    /// # Summary
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is true. If false, nothing happens to the `target`.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIf
    operation ApplyIfC<'T> (op : ('T => Unit : Controlled), bit : Bool, target : 'T) : Unit
    {
        body (...)
        {
            if (bit)
            {
                op(target);
            }
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is true. If false, nothing happens to the `target`.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIf
    operation ApplyIfA<'T> (op : ('T => Unit : Adjoint), bit : Bool, target : 'T) : Unit
    {
        body (...)
        {
            if (bit)
            {
                op(target);
            }
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is true. If false, nothing happens to the `target`.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIf
    operation ApplyIfCA<'T> (op : ('T => Unit : Controlled, Adjoint), bit : Bool, target : 'T) : Unit
    {
        body (...)
        {
            if (bit)
            {
                op(target);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Given an operation op, returns a new operation which
    /// applies the op if a classical control bit is true. If false, nothing happens.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    ///
    /// # Output
    /// A new operation which is op if the classical control bit is true.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CControlledC
    /// - Microsoft.Quantum.Canon.CControlledA
    /// - Microsoft.Quantum.Canon.CControlledCA
    function CControlled<'T> (op : ('T => Unit)) : ((Bool, 'T) => Unit)
    {
        return ApplyIf(op, _, _);
    }
    
    
    /// # Summary
    /// Given an operation op, returns a new operation which
    /// applies the op if a classical control bit is true. If false, nothing happens.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    ///
    /// # Output
    /// A new operation which is op if the classical control bit is true.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CControlled
    function CControlledC<'T> (op : ('T => Unit : Controlled)) : ((Bool, 'T) => Unit : Controlled)
    {
        return ApplyIfC(op, _, _);
    }
    
    
    /// # Summary
    /// Given an operation op, returns a new operation which
    /// applies the op if a classical control bit is true. If false, nothing happens.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    ///
    /// # Output
    /// A new operation which is op if the classical control bit is true.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CControlled
    function CControlledA<'T> (op : ('T => Unit : Adjoint)) : ((Bool, 'T) => Unit : Adjoint)
    {
        return ApplyIfA(op, _, _);
    }
    
    
    /// # Summary
    /// Given an operation op, returns a new operation which
    /// applies the op if a classical control bit is true. If false, nothing happens.
    /// The modifier `CA` indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    ///
    /// # Output
    /// A new operation which is op if the classical control bit is true.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CControlled
    function CControlledCA<'T> (op : ('T => Unit : Controlled, Adjoint)) : ((Bool, 'T) => Unit : Controlled, Adjoint)
    {
        return ApplyIfCA(op, _, _);
    }
    
}


