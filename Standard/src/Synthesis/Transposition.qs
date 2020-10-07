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
    /// Swaps two amplitudes in a state vector.
    ///
    /// # Description
    /// This operation swaps the amplitude at index `a` with the
    /// amplitude at index `b` in the given state-vector of
    /// `register` of length $n$.  If `a` equals `b`, the state-vector
    /// is not changed.
    ///
    /// # Input
    /// ## a
    /// First index (must be a value from 0 to $2^n - 1$)
    /// ## b
    /// Second index (must be a value from 0 to $2^n - 1$)
    /// ## qubits
    /// A list of $n$ qubits to which the transposition is applied to.
    ///
    /// # Example
    /// Prepare a uniform superposition of number states $|1\rangle$, $|2\rangle$, and
    /// $|3\rangle$ on 2 qubits.
    /// ```Q#
    /// using (qubits = Qubit[2]) {
    ///   let register = LittleEndian(qubits);
    ///   PrepareUniformSuperposition(3, register);
    ///   ApplyTransposition(0, 3, register);
    /// }
    /// ```
    operation ApplyTransposition(a : Int, b : Int, qubits : LittleEndian) : Unit is Adj + Ctl {
        let qs = qubits!;
        let n = Length(qs);

        Fact(a >= 0 and a < 2^n, $"Argument a must be value from 0 to {2^n - 1}");
        Fact(b >= 0 and b < 2^n, $"Argument b must be value from 0 to {2^n - 1}");

        if (a != b) {
            let abits = IntAsBoolArray(a, n);
            let bbits = IntAsBoolArray(b, n);
            let diff = IntegerBits(a ^^^ b, n);

            let BitControlledX = ControlledOnBitString(_, X);

            within {
                for (target in Most(diff)) {
                    (BitControlledX(bbits[...target - 1] + abits[target + 1...]))(Excluding([target], qs), qs[target]);
                }
            } apply {
                let target = Tail(diff);
                (BitControlledX(bbits[...target - 1] + abits[target + 1...]))(Excluding([target], qs), qs[target]);
            }
        }
    }
}
