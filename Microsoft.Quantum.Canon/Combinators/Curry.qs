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

    // TODO: functor variants of the above

}