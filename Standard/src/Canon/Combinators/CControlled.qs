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
    operation ApplyIf<'T> (op : ('T => Unit), bit : Bool, target : 'T) : Unit {
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
    operation ApplyIfC<'T> (op : ('T => Unit is Ctl), bit : Bool, target : 'T) : Unit {
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
    operation ApplyIfA<'T> (op : ('T => Unit is Adj), bit : Bool, target : 'T) : Unit {
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
    operation ApplyIfCA<'T> (op : ('T => Unit is Ctl + Adj), bit : Bool, target : 'T) : Unit {
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
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is One. If the Result value is Zero, nothing happens.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
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
    operation ApplyIfOne<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit {
        if (result == One)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is One. If the Result value is Zero, nothing happens.
    /// The modifier C indicates that the operation is controllable.
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
    operation ApplyIfOneC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl {
        if (result == One)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is One. If the Result value is Zero, nothing happens.
    /// The modifier A indicates that the operation is adjointable.
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
    operation ApplyIfOneA<'T> (result : Result, (op : ('T => Unit is Adj + Ctl), target : 'T)) : Unit is Adj + Ctl {
        if (result == One)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is One. If the Result value is Zero, nothing happens.
    /// The modifier CA indicates that the operation is controllable and adjointable.
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
    operation ApplyIfOneCA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj {
        if (result == One)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is Zero. If the Result value is One, nothing happens.
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
    operation ApplyIfZero<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit {
        if (result == Zero)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is Zero. If the Result value is One, nothing happens.
    /// The modifier C indicates that the operation is controllable.
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
    operation ApplyIfZeroC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl {
        if (result == Zero)
        {
            op(target);
        }
    }
    
    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is Zero. If the Result value is One, nothing happens.
    /// The modifier A indicates that the operation is adjointable.
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
    operation ApplyIfZeroA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj {
        if (result == Zero)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given an operation and a Result value, applies the operation
    /// if the Result value is Zero. If the Result value is One, nothing happens.
    /// The modifier CA indicates that the operation is controllable and adjointable.
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
    operation ApplyIfZeroCA<'T> (result : Result, (op : ('T => Unit is Ctl + Adj), target : 'T)) : Unit is Ctl + Adj {
        if (result == Zero)
        {
            op(target);
        }
    }

    /// # Summary
    /// Given two operations and a Result value, applies the first operation
    /// if the Result value is Zero. If the Result value is One, the second operation is applied.
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
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseR<'T,'U> (result : Result, (opZero : ('T => Unit), targetZero : 'T), (opOne : ('U => Unit), targetOne : 'U)) : Unit {
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
    /// Given two operations and a Result value, applies the first operation
    /// if the Result value is Zero. If the Result value is One, the second operation is applied.
    /// The modifier A indicates that the operation is adjointable.
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
    /// - Microsoft.Quantum.Canon.ApplyIfElseR
    operation ApplyIfElseRA<'T,'U> (result : Result, (opZero : ('T => Unit is Adj), targetZero : 'T), (opOne : ('U => Unit is Adj), targetOne : 'U)) : Unit is Adj {
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
    /// Given two operations and a Result value, applies the first operation
    /// if the Result value is Zero. If the Result value is One, the second operation is applied.
    /// The modifier C indicates that the operation is controllable.
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
    /// - Microsoft.Quantum.Canon.ApplyIfElseR
    operation ApplyIfElseRC<'T,'U> (result : Result, (opZero : ('T => Unit is Ctl), targetZero : 'T), (opOne : ('U => Unit is Ctl), targetOne : 'U)) : Unit is Ctl {
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
    /// Given two operations and a Result value, applies the first operation
    /// if the Result value is Zero. If the Result value is One, the second operation is applied.
    /// The modifier CA indicates that the operation is controllable and adjointable.
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
    /// - Microsoft.Quantum.Canon.ApplyIfElseR
    operation ApplyIfElseRCA<'T,'U> (result : Result, (opZero : ('T => Unit is Ctl + Adj), targetZero : 'T), (opOne : ('U => Unit is Ctl + Adj), targetOne : 'U)) : Unit is Ctl + Adj {
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
    /// Given two operations and a Boolean value, applies the first operation
    /// if the Boolean value is True. If the Boolean value is False, the second operation is applied.
    ///
    /// # Input
    /// ## bit
    /// a measurement bit that controls whether op is applied or not.
    /// ## opTrue
    /// An operation to be conditionally applied if bit is True
    /// ## targetTrue
    /// The input to which the operation `opTrue` is applied.
    /// ## opFalse
    /// An operation to be conditionally applied if bit is False
    /// ## targetFalse
    /// The input to which the operation `opFalse` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opTrue` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `opFalse` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfTrue
    /// - Microsoft.Quantum.Canon.ApplyIfFalse
    /// - Microsoft.Quantum.Canon.ApplyIfElseBC
    /// - Microsoft.Quantum.Canon.ApplyIfElseBA
    /// - Microsoft.Quantum.Canon.ApplyIfElseBCA
    operation ApplyIfElseB<'T,'U> (bit : Bool, (opTrue : ('T => Unit), targetTrue : 'T), (opFalse : ('U => Unit), targetFalse : 'U)) : Unit {
        if (bit)
        {
            opTrue(targetTrue);
        }
        else
        {
            opFalse(targetFalse);
        }
    }

    /// # Summary
    /// Given two operations and a Boolean value, applies the first operation
    /// if the Boolean value is True. If the Boolean value is False, the second operation is applied.
    /// The modifier A indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## bit
    /// a measurement bit that controls whether op is applied or not.
    /// ## opTrue
    /// An operation to be conditionally applied if bit is True
    /// ## targetTrue
    /// The input to which the operation `opTrue` is applied.
    /// ## opFalse
    /// An operation to be conditionally applied if bit is False
    /// ## targetFalse
    /// The input to which the operation `opFalse` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opTrue` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `opFalse` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElseB
    operation ApplyIfElseBA<'T,'U> (bit : Bool, (opTrue : ('T => Unit is Adj), targetTrue : 'T), (opFalse : ('U => Unit is Adj), targetFalse : 'U)) : Unit is Adj {
        if (bit)
        {
            opTrue(targetTrue);
        }
        else
        {
            opFalse(targetFalse);
        }
    }

    /// # Summary
    /// Given two operations and a Boolean value, applies the first operation
    /// if the Boolean value is True. If the Boolean value is False, the second operation is applied.
    /// The modifier C indicates that the operation is controllable.
    ///
    /// # Input
    /// ## bit
    /// a measurement bit that controls whether op is applied or not.
    /// ## opTrue
    /// An operation to be conditionally applied if bit is True
    /// ## targetTrue
    /// The input to which the operation `opTrue` is applied.
    /// ## opFalse
    /// An operation to be conditionally applied if bit is False
    /// ## targetFalse
    /// The input to which the operation `opFalse` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opTrue` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `opFalse` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElseB
    operation ApplyIfElseBC<'T,'U> (bit : Bool, (opTrue : ('T => Unit is Ctl), targetTrue : 'T), (opFalse : ('U => Unit is Ctl), targetFalse : 'U)) : Unit is Ctl {
        if (bit)
        {
            opTrue(targetTrue);
        }
        else
        {
            opFalse(targetFalse);
        }
    }

    /// # Summary
    /// Given two operations and a Boolean value, applies the first operation
    /// if the Boolean value is True. If the Boolean value is False, the second operation is applied.
    /// The modifier CA indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## bit
    /// a measurement bit that controls whether op is applied or not.
    /// ## opTrue
    /// An operation to be conditionally applied if bit is True
    /// ## targetTrue
    /// The input to which the operation `opTrue` is applied.
    /// ## opFalse
    /// An operation to be conditionally applied if bit is False
    /// ## targetFalse
    /// The input to which the operation `opFalse` is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `opTrue` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `opFalse` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfElseB
    operation ApplyIfElseBCA<'T,'U> (bit : Bool, (opTrue : ('T => Unit is Ctl + Adj), targetTrue : 'T), (opFalse : ('U => Unit is Ctl + Adj), targetFalse : 'U)) : Unit is Ctl + Adj {
        if (bit)
        {
            opTrue(targetTrue);
        }
        else
        {
            opFalse(targetFalse);
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
    function CControlled<'T> (op : ('T => Unit)) : ((Bool, 'T) => Unit) {
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
    function CControlledC<'T> (op : ('T => Unit is Ctl)) : ((Bool, 'T) => Unit is Ctl) {
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
    function CControlledA<'T> (op : ('T => Unit is Adj)) : ((Bool, 'T) => Unit is Adj) {
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
    function CControlledCA<'T> (op : ('T => Unit is Ctl + Adj)) : ((Bool, 'T) => Unit is Ctl + Adj) {
        return ApplyIfCA(op, _, _);
    }
    
}


