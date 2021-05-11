// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

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
    /// Argument used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    /// ## 'U
    /// The return type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.DelayC
    /// - Microsoft.Quantum.Canon.DelayA
    /// - Microsoft.Quantum.Canon.DelayCA
    /// - Microsoft.Quantum.Canon.Delayed
    operation Delay<'T, 'U> (op : ('T => 'U), arg : 'T, aux : Unit) : 'U {
        return op(arg);
    }


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
    /// Argument used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay
    /// - Microsoft.Quantum.Canon.Delayed
    operation DelayA<'T> ( op : ('T => Unit is Adj), arg : 'T, aux : Unit) : Unit is Adj {
        op(arg);
    }


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
    /// Argument used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay 
    /// - Microsoft.Quantum.Canon.Delayed
    operation DelayC<'T> ( op : ('T => Unit is Ctl), arg : 'T, aux : Unit) : Unit is Ctl {
        op(arg);
    }


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
    /// Argument used to delay the application of operation by using 
    /// partial application.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Delay
    /// - Microsoft.Quantum.Canon.Delayed
    operation DelayCA<'T> ( op : ('T => Unit is Ctl + Adj), arg : 'T, aux : Unit) : Unit is Ctl + Adj {
        op(arg);
    }


    /// # Summary
    /// Returns an operation that applies 
    /// given operation with given argument.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## arg
    /// The input to which the operation is applied.
    ///
    /// # Output
    /// A new operation which applies `op` with input `arg`
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be delayed.
    /// ## 'U
    /// The return type of the operation to be delayed.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.DelayedC
    /// - Microsoft.Quantum.Canon.DelayedA
    /// - Microsoft.Quantum.Canon.DelayedCA
    /// - Microsoft.Quantum.Canon.Delay
    function Delayed<'T, 'U> ( op : ('T => 'U), arg : 'T) : (Unit => 'U) {
        return Delay(op, arg, _);
    }


    /// # Summary
    /// Returns an operation that applies 
    /// given operation with given argument.
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
    /// - Microsoft.Quantum.Canon.Delay
    function DelayedA<'T> ( op : ('T => Unit is Adj), arg : 'T) : (Unit => Unit is Adj) {
        return DelayA(op, arg, _);
    }


    /// # Summary
    /// Returns an operation that applies 
    /// given operation with given argument.
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
    /// - Microsoft.Quantum.Canon.Delay
    function DelayedC<'T> ( op : ('T => Unit is Ctl), arg : 'T) : (Unit => Unit is Ctl) {
        return DelayC(op, arg, _);
    }


    /// # Summary
    /// Returns an operation that applies 
    /// given operation with given argument.
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
    /// - Microsoft.Quantum.Canon.Delay
    function DelayedCA<'T> ( op : ('T => Unit is Ctl + Adj), arg : 'T) : (Unit => Unit is Ctl + Adj) {
        return DelayCA(op, arg, _);
    }
}
