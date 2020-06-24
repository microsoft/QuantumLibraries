// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    internal operation ApplyBound<'T> (operations : ('T => Unit)[], target : 'T) : Unit {
        for (op in operations) {
            op(target);
        }
    }

    /// # Summary
    /// Given an array of operations acting on a single input,
    /// produces a new operation that
    /// performs each given operation in sequence.
    ///
    /// # Input
    /// ## operations
    /// A sequence of operations to be performed on a given input.
    ///
    /// # Output
    /// A new operation that performs each given operation in sequence
    /// on its input.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations in the array act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// let bound = Bound([U, V]);
    /// bound(x);
    /// ```
    /// and
    /// ```qsharp
    /// U(x); V(x);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundC
    /// - Microsoft.Quantum.Canon.BoundA
    /// - Microsoft.Quantum.Canon.BoundCA
    function Bound<'T> (operations : ('T => Unit)[]) : ('T => Unit) {
        return ApplyBound(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundA
    internal operation ApplyBoundA<'T> (operations : ('T => Unit is Adj)[], target : 'T)
    : Unit is Adj {
        for (op in operations) {
            op(target);
        }
    }


    /// # Summary
    /// Given an array of operations acting on a single input,
    /// produces a new operation that
    /// performs each given operation in sequence.
    /// The modifier `A` indicates that all operations in the array are adjointable.
    ///
    /// # Input
    /// ## operations
    /// A sequence of operations to be performed on a given input.
    ///
    /// # Output
    /// A new operation that performs each given operation in sequence
    /// on its input.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations in the array act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// let bound = BoundA([U, V]);
    /// bound(x);
    /// ```
    /// and
    /// ```qsharp
    /// U(x); V(x);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    function BoundA<'T> (operations : ('T => Unit is Adj)[])
    : ('T => Unit is Adj) {
        return ApplyBoundA(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundC
    internal operation ApplyBoundC<'T> (operations : ('T => Unit is Ctl)[], target : 'T)
    : Unit is Ctl {
        for (op in operations) {
            op(target);
        }
    }


    /// # Summary
    /// Given an array of operations acting on a single input,
    /// produces a new operation that
    /// performs each given operation in sequence.
    /// The modifier `C` indicates that all operations in the array are controllable.
    ///
    /// # Input
    /// ## operations
    /// A sequence of operations to be performed on a given input.
    ///
    /// # Output
    /// A new operation that performs each given operation in sequence
    /// on its input.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations in the array act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// let bound = BoundC([U, V]);
    /// bound(x);
    /// ```
    /// and
    /// ```qsharp
    /// U(x); V(x);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    function BoundC<'T> (operations : ('T => Unit is Ctl)[]) : ('T => Unit is Ctl) {
        return ApplyBoundC(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundCA
    internal operation ApplyBoundCA<'T> (operations : ('T => Unit is Adj + Ctl)[], target : 'T)
    : Unit is Adj + Ctl {
        for (op in operations) {
            op(target);
        }
    }


    /// # Summary
    /// Given an array of operations acting on a single input,
    /// produces a new operation that
    /// performs each given operation in sequence.
    /// The modifier `CA` indicates that all operations in the array are adjointable
    /// and controllable.
    ///
    /// # Input
    /// ## operations
    /// A sequence of operations to be performed on a given input.
    ///
    /// # Output
    /// A new operation that performs each given operation in sequence
    /// on its input.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the operations in the array act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// let bound = BoundCA([U, V]);
    /// bound(x);
    /// ```
    /// and
    /// ```qsharp
    /// U(x); V(x);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    function BoundCA<'T> (operations : ('T => Unit is Adj + Ctl)[]) : ('T => Unit is Adj + Ctl) {
        return ApplyBoundCA(operations, _);
    }

}


