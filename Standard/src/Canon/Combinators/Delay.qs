// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{

    /// # Summary
    /// Applies a given operation with a delay.
    ///
    /// # Description
    /// Given an operation and an input to that operation, applies
    /// the operation once an additional input is provided.
    /// In particular, the expression `Delay(op, arg, _)` is an operation that
    /// applies `op` to `arg` when called.
    /// Expression `Delay(op,arg,_)` allows to delay the application of `op`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## arg
    /// The input to which the operation is applied.
    /// ## aux
    /// Argument of type Unit used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.DelayC
    /// - Microsoft.Quantum.Canon.DelayA
    /// - Microsoft.Quantum.Canon.DelayCA    
    operation Delay<'T> ( op : ('T => Unit), arg : 'T, aux : Unit) : Unit
    {
        op(arg);
    }


    /// # Summary
    /// Applies operation `op` with argument `arg`.
    /// Expression `Delay(op,arg,_)` allows to delay the application of `op`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## arg
    /// The input to which the operation is applied.
    /// ## aux
    /// Argument of type Unit used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay    
    operation DelayA<'T> ( op : ('T => Unit is Adj), arg : 'T, aux : Unit) : Unit is Adj
    {
        op(arg);
    }


    /// # Summary
    /// Applies operation `op` with argument `arg`.
    /// Expression `Delay(op,arg,_)` allows to delay the application of `op`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## arg
    /// The input to which the operation is applied.
    /// ## aux
    /// Argument of type Unit used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay 
    operation DelayC<'T> ( op : ('T => Unit is Ctl), arg : 'T, aux : Unit) : Unit is Ctl
    {
        op(arg);
    }


    /// # Summary
    /// Applies operation `op` with argument `arg`.
    /// Expression `Delay(op,arg,_)` allows to delay the application of `op`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## arg
    /// The input to which the operation is applied.
    /// ## aux
    /// Argument of type Unit used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay 
    operation DelayCA<'T> ( op : ('T => Unit is Ctl + Adj), arg : 'T, aux : Unit) : Unit is Ctl + Adj
    {
        op(arg);
    }


    /// # Summary
    /// Returns an operation of type `Unit => Unit` that applies 
    /// operation `op` with argument `arg`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied as a result of applying return value
    /// ## arg
    /// The input to which the operation `op` is applied.
    ///
    /// # Output
    /// A new operation which applies `op` with input `arg`
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.DelayedC
    /// - Microsoft.Quantum.Canon.DelayedA
    /// - Microsoft.Quantum.Canon.DelayedCA    
    function Delayed<'T> ( op : ('T => Unit), arg : 'T) : (Unit => Unit)
    {
        return Delay(op, arg, _);
    }


    /// # Summary
    /// Returns an operation of type `Unit => Unit is Adj` that applies 
    /// operation `op` with argument `arg`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied as a result of applying return value
    /// ## arg
    /// The input to which the operation `op` is applied.
    ///
    /// # Output
    /// A new operation which applies `op` with input `arg`
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delayed 
    function DelayedA<'T> ( op : ('T => Unit is Adj), arg : 'T) : (Unit => Unit is Adj)
    {
        return DelayA(op, arg, _);
    }


    /// # Summary
    /// Returns an operation of type `Unit => Unit is Ctl` that applies 
    /// operation `op` with argument `arg`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied as a result of applying return value
    /// ## arg
    /// The input to which the operation `op` is applied.
    ///
    /// # Output
    /// A new operation which applies `op` with input `arg`
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delayed  
    function DelayedC<'T> ( op : ('T => Unit is Ctl), arg : 'T) : (Unit => Unit is Ctl)
    {
        return DelayC(op, arg, _);
    }


    /// # Summary
    /// Returns an operation of type `Unit => Unit is Ctl + Adj` that applies 
    /// operation `op` with argument `arg`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied as a result of applying return value
    /// ## arg
    /// The input to which the operation `op` is applied.
    ///
    /// # Output
    /// A new operation which applies `op` with input `arg`
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delayed
    function DelayedCA<'T> ( op : ('T => Unit is Ctl + Adj), arg : 'T) : (Unit => Unit is Ctl + Adj)
    {
        return DelayCA(op, arg, _);
    }
}
