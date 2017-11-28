// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    // TODO: remove references to qubits in doc comments.
    // TODO: remove references to 1 in comments (Bind ← Bind1).
    // TODO: update to DFM styled ///.

    operation BindImpl<'T>(operations : ('T => ())[], target : 'T) : () {
        body {
            for (idxOperation in 0..Length(operations) - 1) {
                let op = operations[idxOperation];
                op(target);
            }
        }
    }

    /// # Summary
    /// Given an array of operations acting on a single qubit,
    /// produces a new operation that
    /// performs each given operation in sequence.
    ///
    /// ## Remark
    /// See Bind1A, Bind1C, and Bind1AC functor variants.
    function Bind<'T>(operations : ('T => ())[]) : ('T => ()) {
        return BindImpl(operations, _);
    }

    operation BindAImpl<'T>(operations : ('T => () : Adjoint)[], target : 'T) : () {
        body {
            BindImpl(operations, target);
        }
        adjoint {
            // TODO: replace with an implementation based on Reversed : 'T[] -> 'T[]
            //       and AdjointAll : ('T => () : Adjointable)[] -> ('T => () : Adjointable).
            for (idxOperation in Length(operations) - 1..0) {
                let op = (Adjoint operations[idxOperation]);
                op(target);
            }
        }
    }

    /// # Summary
    /// Given an array of operations acting on a single qubit,
    /// produces a new operation that
    /// performs each given operation in sequence.
    ///
    /// # Remarks
    /// See Bind1, Bind1C, and Bind1AC functor variants.
    function BindA<'T>(operations : ('T => () : Adjoint)[]) : ('T => () : Adjoint) {
        return BindAImpl(operations, _);
    }

    operation BindCImpl<'T>(operations : ('T => () : Controlled)[], target : 'T) : () {
        body {
            BindImpl(operations, target);
        }

        controlled (controls) {
            for (idxOperation in 0..Length(operations) - 1) {
                let op = (Controlled operations[idxOperation]);
                op(controls, target);
            }
        }
    }

    /// # Summary
    /// Given an array of operations acting on a single qubit,
    /// produces a new operation that
    /// performs each given operation in sequence.
    ///
    /// # Remarks
    /// See Bind1, Bind1A, and Bind1AC functor variants.
    function BindC<'T>(operations : ('T => () : Controlled)[]) : ('T => () : Controlled) {
        return BindCImpl(operations, _);
    }

    operation BindACImpl<'T>(operations : ('T => () : Adjoint, Controlled)[], target : 'T) : () {
        body {
            BindImpl(operations, target);
        }

        adjoint {
            (Adjoint BindAImpl)(operations, target);
        }
        controlled (controls) {
            (Controlled BindCImpl)(controls, (operations, target));
        }

        controlled adjoint (controls) {
            for (idxOperation in Length(operations) - 1..0) {
                let op = (Controlled Adjoint operations[idxOperation]);
                op(controls, target);
            }
        }
    }

    /// # Summary
    /// Given an array of operations acting on a single qubit,
    /// produces a new operation that
    /// performs each given operation in sequence.
    ///
    /// # Remarks
    /// See Bind1, Bind1A, and Bind1AC functor variants.
    function BindAC<'T>(operations : ('T => () : Adjoint, Controlled)[]) : ('T => () : Adjoint, Controlled) {
        return BindACImpl(operations, _);
    }

}
