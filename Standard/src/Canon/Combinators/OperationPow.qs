// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    internal operation ApplyOperationRepeatedly<'T> (op : ('T => Unit), power : Int, target : 'T)
    : Unit {
        for (idxApplication in 0 .. power - 1) {
            op(target);
        }
    }


    internal operation ApplyOperationRepeatedlyC<'T> (op : ('T => Unit is Ctl), power : Int, target : 'T)
    : Unit is Ctl {
        for (idxApplication in 0 .. power - 1) {
            op(target);
        }
    }


    internal operation ApplyOperationRepeatedlyA<'T> (op : ('T => Unit is Adj), power : Int, target : 'T)
    : Unit is Adj {
        for (idxApplication in 0 .. power - 1) {
            op(target);
        }
    }


    internal operation ApplyOperationRepeatedlyCA<'T> (op : ('T => Unit is Adj + Ctl), power : Int, target : 'T)
    : Unit is Adj + Ctl {
        for (idxApplication in 0 .. power - 1) {
            op(target);
        }
    }


    /// # Summary
    /// Raises an operation to a power.
    ///
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## op
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.OperationPowC
    /// - Microsoft.Quantum.Canon.OperationPowA
    /// - Microsoft.Quantum.Canon.OperationPowCA
    function OperationPow<'T> (op : ('T => Unit), power : Int) : ('T => Unit) {
        return ApplyOperationRepeatedly(op, power, _);
    }


    /// # Summary
    /// Raises an operation to a power.
    /// The modifier `C` indicates that the operation is controllable.
    ///
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## op
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.OperationPow
    function OperationPowC<'T> (op : ('T => Unit is Ctl), power : Int) : ('T => Unit is Ctl) {
        return ApplyOperationRepeatedlyC(op, power, _);
    }


    /// # Summary
    /// Raises an operation to a power.
    /// The modifier `A` indicates that the operation is adjointable.
    ///
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## op
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.OperationPow
    function OperationPowA<'T> (op : ('T => Unit is Adj), power : Int) : ('T => Unit is Adj) {
        return ApplyOperationRepeatedlyA(op, power, _);
    }


    /// # Summary
    /// Raises an operation to a power.
    /// The modifier `A` indicates that the operation is controllable and adjointable.
    ///
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## op
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.OperationPow
    function OperationPowCA<'T> (op : ('T => Unit is Ctl + Adj), power : Int) : ('T => Unit is Ctl + Adj) {
        return ApplyOperationRepeatedlyCA(op, power, _);
    }

}
