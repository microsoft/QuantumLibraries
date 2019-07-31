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
    operation ApplyIfC<'T> (op : ('T => Unit is Ctl), bit : Bool, target : 'T) : Unit
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
    operation ApplyIfA<'T> (op : ('T => Unit is Adj), bit : Bool, target : 'T) : Unit
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
    operation ApplyIfCA<'T> (op : ('T => Unit is Ctl + Adj), bit : Bool, target : 'T) : Unit
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
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is One. If result is Zero, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOneC
    /// - Microsoft.Quantum.Canon.ApplyIfOneA
    /// - Microsoft.Quantum.Canon.ApplyIfOneCA
    operation ApplyIfOne<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit
    {
        if (result == One)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is One. If result is Zero, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    operation ApplyIfOneC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl
    {
        if (result == One)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is One. If result is Zero, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    operation ApplyIfOneA<'T> (result : Result, (op : ('T => Unit is Adj + Ctl), target : 'T)) : Unit is Adj + Ctl
    {
        if (result == One)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is One. If result is Zero, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    operation ApplyIfOneCA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj
    {
        if (result == One)
        {
            op(target);
        }
    }

	/// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is Zero. If result is One, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZeroC
    /// - Microsoft.Quantum.Canon.ApplyIfZeroA
    /// - Microsoft.Quantum.Canon.ApplyIfZeroCA
    operation ApplyIfZero<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit
    {
        if (result == Zero)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is Zero. If result is One, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    operation ApplyIfZeroC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl
    {
        if (result == Zero)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is Zero. If result is One, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    operation ApplyIfZeroA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj
    {
        if (result == Zero)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given an operation `op` and a Result value `result`, applies `op` to the `target`
    /// if `result` is Zero. If result is One, nothing happens to the `target`.
    ///
    /// # Input
	/// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## op
    /// An operation to be conditionally applied.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    operation ApplyIfZeroCA<'T> (result : Result, (op : ('T => Unit is Ctl + Adj), target : 'T)) : Unit is Ctl + Adj
    {
        if (result == Zero)
        {
            op(target);
        }
    }

	/// # Summary
    /// Given operations `opZero`, `opOne` and a Result value `result`, applies `op` to the `targetZero`
    /// if `result` is Zero. If result is One, applies `opOne` to the `targetOne`.
    ///
    /// # Input
	/// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## opZero
    /// An operation to be conditionally applied if result is Zero
    /// ## targetZero
    /// The input to which the operation `opZero` is applied.
	/// ## opOne
    /// An operation to be conditionally applied if result is One
    /// ## targetOne
    /// The input to which the operation `opOne` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opZero` to be conditionally applied.
	/// ## 'U
    /// The input type of the operation `opOne` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElseC
    /// - Microsoft.Quantum.Canon.ApplyIfElseA
    /// - Microsoft.Quantum.Canon.ApplyIfElseCA
    operation ApplyIfElse<'T,'U> (result : Result, (opZero : ('T => Unit), targetZero : 'T), (opOne : ('U => Unit), targetOne : 'U)) : Unit
    {
        if (result == Zero)
        {
            opZero(targetZero);
        }
		else
		{
			opOne(targetOne);
		}
    }

	/// # Summary
    /// Given operations `opZero`, `opOne` and a Result value `result`, applies `op` to the `targetZero`
    /// if `result` is Zero. If result is One, applies `opOne` to the `targetOne`.
    ///
    /// # Input
	/// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## opZero
    /// An operation to be conditionally applied if result is Zero
    /// ## targetZero
    /// The input to which the operation `opZero` is applied.
	/// ## opOne
    /// An operation to be conditionally applied if result is One
    /// ## targetOne
    /// The input to which the operation `opOne` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opZero` to be conditionally applied.
	/// ## 'U
    /// The input type of the operation `opOne` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElse
    operation ApplyIfElseA<'T,'U> (result : Result, (opZero : ('T => Unit is Adj), targetZero : 'T), (opOne : ('U => Unit is Adj), targetOne : 'U)) : Unit is Adj
    {
        if (result == Zero)
        {
            opZero(targetZero);
        }
		else
		{
			opOne(targetOne);
		}
    }

	/// # Summary
    /// Given operations `opZero`, `opOne` and a Result value `result`, applies `op` to the `targetZero`
    /// if `result` is Zero. If result is One, applies `opOne` to the `targetOne`.
    ///
    /// # Input
	/// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## opZero
    /// An operation to be conditionally applied if result is Zero
    /// ## targetZero
    /// The input to which the operation `opZero` is applied.
	/// ## opOne
    /// An operation to be conditionally applied if result is One
    /// ## targetOne
    /// The input to which the operation `opOne` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opZero` to be conditionally applied.
	/// ## 'U
    /// The input type of the operation `opOne` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElse
    operation ApplyIfElseC<'T,'U> (result : Result, (opZero : ('T => Unit is Ctl), targetZero : 'T), (opOne : ('U => Unit is Ctl), targetOne : 'U)) : Unit is Ctl
    {
        if (result == Zero)
        {
            opZero(targetZero);
        }
		else
		{
			opOne(targetOne);
		}
    }

	/// # Summary
    /// Given operations `opZero`, `opOne` and a Result value `result`, applies `op` to the `targetZero`
    /// if `result` is Zero. If result is One, applies `opOne` to the `targetOne`.
    ///
    /// # Input
	/// ## result
    /// a measurement result that controls whether op is applied or not.
    /// ## opZero
    /// An operation to be conditionally applied if result is Zero
    /// ## targetZero
    /// The input to which the operation `opZero` is applied.
	/// ## opOne
    /// An operation to be conditionally applied if result is One
    /// ## targetOne
    /// The input to which the operation `opOne` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opZero` to be conditionally applied.
	/// ## 'U
    /// The input type of the operation `opOne` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElse
    operation ApplyIfElseCA<'T,'U> (result : Result, (opZero : ('T => Unit is Ctl + Adj), targetZero : 'T), (opOne : ('U => Unit is Ctl + Adj), targetOne : 'U)) : Unit is Ctl + Adj
    {
        if (result == Zero)
        {
            opZero(targetZero);
        }
		else
		{
			opOne(targetOne);
		}
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
    function CControlledC<'T> (op : ('T => Unit is Ctl)) : ((Bool, 'T) => Unit is Ctl)
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
    function CControlledA<'T> (op : ('T => Unit is Adj)) : ((Bool, 'T) => Unit is Adj)
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
    function CControlledCA<'T> (op : ('T => Unit is Ctl + Adj)) : ((Bool, 'T) => Unit is Ctl + Adj)
    {
        return ApplyIfCA(op, _, _);
    }
    
}


