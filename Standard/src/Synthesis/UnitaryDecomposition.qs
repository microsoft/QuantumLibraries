// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;


    // # Summary
    /// Applies single-qubit gate defined by 2x2 unitary matrix.
    ///
    /// # Description
    /// This operation swaps the amplitude at index `a` with the
    /// amplitude at index `b` in the given state-vector of
    /// `register` of length $n$.  If `a` equals `b`, the state-vector
    /// is not changed.
    ///
    /// # Input
    /// ## unitary
    /// 2x2 unitary matrix describing the operation.
    /// ## qubit
    /// Qubit to which apply the operation.
    operation ApplySingleQubitUnitary(qubit : Qubit) : Unit is Adj + Ctl {
        X(qubit);
    }


    operation ApplyUnitary(qubits : LittleEndian) : Unit is Adj + Ctl {
        ApplySingleQubitUnitary(qubits[0]);
    }    
}