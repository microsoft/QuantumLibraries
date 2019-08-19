// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Given an operation op, returns a new operation which
    /// applies the op if a classical control bit is true. If `false`, nothing happens.
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
    /// applies the op if a classical control bit is true. If `false`, nothing happens.
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
    /// applies the op if a classical control bit is true. If `false`, nothing happens.
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
    /// applies the op if a classical control bit is true. If `false`, nothing happens.
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

