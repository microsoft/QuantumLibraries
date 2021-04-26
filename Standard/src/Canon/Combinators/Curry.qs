// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    internal function WithFirstInputApplied<'T, 'U> (op : (('T, 'U) => Unit), arg1 : 'T) : ('U => Unit) {
        return op(arg1, _);
    }

    internal function WithFirstInputAppliedA<'T, 'U> (op : (('T, 'U) => Unit is Adj), arg1 : 'T) : ('U => Unit is Adj) {
        return op(arg1, _);
    }

    internal function WithFirstInputAppliedC<'T, 'U> (op : (('T, 'U) => Unit is Ctl), arg1 : 'T) : ('U => Unit is Ctl) {
        return op(arg1, _);
    }

    internal function WithFirstInputAppliedCA<'T, 'U> (op : (('T, 'U) => Unit is Adj + Ctl), arg1 : 'T) : ('U => Unit is Adj + Ctl) {
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CurriedOpC
    /// - Microsoft.Quantum.Canon.CurriedOpA
    /// - Microsoft.Quantum.Canon.CurriedOpCA
    function CurriedOp<'T, 'U> (op : (('T, 'U) => Unit)) : ('T -> ('U => Unit)) {
        return WithFirstInputApplied(op, _);
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CurriedOp
    /// - Microsoft.Quantum.Canon.CurriedOpC
    /// - Microsoft.Quantum.Canon.CurriedOpCA
    function CurriedOpA<'T, 'U> (op : (('T, 'U) => Unit is Adj)) : ('T -> ('U => Unit is Adj)) {
        return WithFirstInputAppliedA(op, _);
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CurriedOp
    /// - Microsoft.Quantum.Canon.CurriedOpA
    /// - Microsoft.Quantum.Canon.CurriedOpCA
    function CurriedOpC<'T, 'U> (op : (('T, 'U) => Unit is Ctl)) : ('T -> ('U => Unit is Ctl)) {
        return WithFirstInputAppliedC(op, _);
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CurriedOp
    /// - Microsoft.Quantum.Canon.CurriedOpC
    /// - Microsoft.Quantum.Canon.CurriedOpA
    function CurriedOpCA<'T, 'U> (op : (('T, 'U) => Unit is Adj + Ctl)) : ('T -> ('U => Unit is Adj + Ctl)) {
        return WithFirstInputAppliedCA(op, _);
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
    /// - Microsoft.Quantum.Canon.UncurriedOpC
    /// - Microsoft.Quantum.Canon.UncurriedOpA
    /// - Microsoft.Quantum.Canon.UncurriedOpCA
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
    /// - Microsoft.Quantum.Canon.UncurriedOp
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
    /// - Microsoft.Quantum.Canon.UncurriedOp
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
    /// - Microsoft.Quantum.Canon.UncurriedOp
    function UncurriedOpCA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl + Adj))) : (('T, 'U) => Unit is Ctl + Adj) {
        return ApplyCurriedOpCA(curriedOp, _, _);
    }

}


