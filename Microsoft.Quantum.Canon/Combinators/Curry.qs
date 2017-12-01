// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    function CurryOpImpl<'T, 'U>(op : (('T, 'U) => ()), arg1 : 'T) : ('U => ()) {
        return op(arg1, _);
    }

    /// # Summary
    /// Given an operation with two inputs, applies Curry's isomorphism
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
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// op(x, y);
    ///
    /// let curried = CurryOp(op);
    /// let partial = curried(x);
    /// partial(y);
    /// ```
    function CurryOp<'T, 'U>(op : (('T, 'U) => ())) : ('T -> ('U => ())) {
        return CurryOpImpl(op, _);
    }

    operation UncurryOpImpl<'T, 'U>(curriedOp : ('T -> ('U => ())), first : 'T, second : 'U) : () {
        body {
            let innerOp = curriedOp(first);
            innerOp(second);
        }
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
    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryopc"
    /// - @"microsoft.quantum.canon.uncurryopa"
    /// - @"microsoft.quantum.canon.uncurryopca"
    function UncurryOp<'T, 'U>(curriedOp : ('T -> ('U => ()))) : (('T, 'U) => ()) {
        return UncurryOpImpl(curriedOp, _, _);
    }

    operation UncurryOpCImpl<'T, 'U>(curriedOp : ('T -> ('U => () : Controlled)), first : 'T, second : 'U) : () {
        body {
            let innerOp = curriedOp(first);
            innerOp(second);
        }

        controlled auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurryOpC<'T, 'U>(curriedOp : ('T -> ('U => () : Controlled))) : (('T, 'U) => () : Controlled) {
        return UncurryOpCImpl(curriedOp, _, _);
    }

    operation UncurryOpAImpl<'T, 'U>(curriedOp : ('T -> ('U => () : Adjoint)), first : 'T, second : 'U) : () {
        body {
            let innerOp = curriedOp(first);
            innerOp(second);
        }

        adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurryOpA<'T, 'U>(curriedOp : ('T -> ('U => () : Adjoint))) : (('T, 'U) => () : Adjoint) {
        return UncurryOpAImpl(curriedOp, _, _);
    }

    operation UncurryOpCAImpl<'T, 'U>(curriedOp : ('T -> ('U => () : Controlled, Adjoint)), first : 'T, second : 'U) : () {
        body {
            let innerOp = curriedOp(first);
            innerOp(second);
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.uncurryop"
    function UncurryOpCA<'T, 'U>(curriedOp : ('T -> ('U => () : Controlled, Adjoint))) : (('T, 'U) => () : Controlled, Adjoint) {
        return UncurryOpCAImpl(curriedOp, _, _);
    }

}
