// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    operation _Bound<'T> (operations : ('T => Unit)[], target : 'T) : Unit {
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
        return _Bound(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundA
    operation _BoundA<'T> (operations : ('T => Unit is Adj)[], target : 'T) : Unit {
        body (...) {
            for (op in operations) {
                op(target);
            }
        }

        adjoint (...) {
            // TODO: replace with an implementation based on Reversed : 'T[] -> 'T[]
            //       and AdjointAll : ('T => () is Adj)[] -> ('T => () is Adj).
            for (op in Reversed(operations)) {
                Adjoint op(target);
            }
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
    function BoundA<'T> (operations : ('T => Unit is Adj)[]) : ('T => Unit is Adj)
    {
        return _BoundA(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundC
    operation _BoundC<'T> (operations : ('T => Unit is Ctl)[], target : 'T) : Unit {
        body (...) {
            for (op in operations) {
                op(target);
            }
        }

        controlled (controls, ...) {
            for (op in operations) {
                Controlled op(controls, target);
            }
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
        return _BoundC(operations, _);
    }


    /// # See Also
    /// - Microsoft.Quantum.Canon.BoundCA
    operation _BoundCA<'T> (operations : ('T => Unit is Adj + Ctl)[], target : 'T) : Unit {
        body (...) {
            for (idxOperation in operations) {
                op(target);
            }
        }

        adjoint (...) {
            for (op in Reversed(operations)) {
                Adjoint op(target);
            }
        }

        controlled (controls, ...) {
            for (op in operations) {
                Controlled op(controls, target);
            }
        }

        controlled adjoint (controls, ...) {
            for (op in Reversed(operations)) {
                Controlled Adjoint op(controls, target);
            }
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
        return _BoundCA(operations, _);
    }

}


