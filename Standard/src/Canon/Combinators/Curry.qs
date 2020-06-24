// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    internal function WithFirstInputApplied<'T, 'U> (op : (('T, 'U) => Unit), arg1 : 'T) : ('U => Unit) {
        return op(arg1, _);
    }


    /// # Summary
    /// Returns a curried version of an operation on two inputs.
    ///
    /// That is, given an operation with two inputs, this function applies Curry's isomorphism
    /// $f(x, y) \equiv f(x)(y)$ to return an operation of one input which
    /// returns an operation of one input.
    ///
    /// # Input
    /// ## op
    /// An operation whose input is a pair.
    ///
    /// # Output
    /// An operation which accepts the first element of a pair and returns
    /// an operation which accepts as its input the second element of the
    /// original operation's input.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first component of a function defined on pairs.
    /// ## 'U
    /// The type of the second component of a function defined on pairs.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// op(x, y);
    ///
    /// let curried = CurriedOp(op);
    /// let partial = curried(x);
    /// partial(y);
    /// ```
    function CurriedOp<'T, 'U> (op : (('T, 'U) => Unit)) : ('T -> ('U => Unit)) {
        return WithFirstInputApplied(op, _);
    }

    internal operation ApplyCurriedOp<'T, 'U> (curriedOp : ('T -> ('U => Unit)), first : 'T, second : 'U) : Unit {
        let innerOp = curriedOp(first);
        innerOp(second);
    }


    /// # Summary
    /// Given a function which returns operations,
    /// returns a new operation which takes both inputs
    /// as a tuple.
    ///
    /// # Input
    /// ## curriedOp
    /// A function which returns operations.
    ///
    /// # Output
    /// A new operation `op` such that `op(t, u)` is equivalent
    /// to `(curriedOp(t))(u)`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first input to a curried operation.
    /// ## 'U
    /// The type of the second input to a curried operation.
    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryopc"
    /// - @"microsoft.quantum.canon.uncurryopa"
    /// - @"microsoft.quantum.canon.uncurryopca"
    function UncurriedOp<'T, 'U> (curriedOp : ('T -> ('U => Unit))) : (('T, 'U) => Unit) {
        return ApplyCurriedOp(curriedOp, _, _);
    }


    internal operation ApplyCurriedOpC<'T, 'U>(curriedOp : ('T -> ('U => Unit is Ctl)), first : 'T, second : 'U)
    : Unit is Ctl {
        let innerOp = curriedOp(first);
        innerOp(second);
    }


    /// # Summary
    /// Given a function which returns operations,
    /// returns a new operation which takes both inputs
    /// as a tuple.
    /// The modifier `C` indicates that the operations are controllable.
    ///
    /// # Input
    /// ## curriedOp
    /// A function which returns operations.
    ///
    /// # Output
    /// A new operation `op` such that `op(t, u)` is equivalent
    /// to `(curriedOp(t))(u)`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first argument of a curried function.
    /// ## 'U
    /// The type of the second argument of a curried function.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurriedOpC<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl)))
    : (('T, 'U) => Unit is Ctl) {
        return ApplyCurriedOpC(curriedOp, _, _);
    }


    internal operation ApplyCurriedOpA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj)), first : 'T, second : 'U)
    : Unit is Adj {
        let innerOp = curriedOp(first);
        innerOp(second);
    }


    /// # Summary
    /// Given a function which returns operations,
    /// returns a new operation which takes both inputs
    /// as a tuple.
    /// The modifier `A` indicates that the operations are adjointable.
    ///
    /// # Input
    /// ## curriedOp
    /// A function which returns operations.
    ///
    /// # Output
    /// A new operation `op` such that `op(t, u)` is equivalent
    /// to `(curriedOp(t))(u)`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first argument of a curried function.
    /// ## 'U
    /// The type of the second argument of a curried function.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurriedOpA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj))) : (('T, 'U) => Unit is Adj) {
        return ApplyCurriedOpA(curriedOp, _, _);
    }


    internal operation ApplyCurriedOpCA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl + Adj)), first : 'T, second : 'U)
    : Unit is Adj + Ctl {
        let innerOp = curriedOp(first);
        innerOp(second);
    }


    /// # Summary
    /// Given a function which returns operations,
    /// returns a new operation which takes both inputs
    /// as a tuple.
    /// The modifier `CA` indicates that the operations are controllable and adjointable.
    ///
    /// # Input
    /// ## curriedOp
    /// A function which returns operations.
    ///
    /// # Output
    /// A new operation `op` such that `op(t, u)` is equivalent
    /// to `(curriedOp(t))(u)`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the first argument of a curried function.
    /// ## 'U
    /// The type of the second argument of a curried function.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurriedOpCA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl + Adj))) : (('T, 'U) => Unit is Ctl + Adj) {
        return ApplyCurriedOpCA(curriedOp, _, _);
    }

}


