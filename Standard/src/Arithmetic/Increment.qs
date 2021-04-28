// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Increments an unsigned quantum register by a classical integer,
    /// using phase rotations.
    ///
    /// # Description
    /// Suppose that `target` encodes an unsigned integer $x$ in a little-endian
    /// encoding and that `increment` is equal to $a$.
    /// Then, this operation implements the unitary $\ket{x} \mapsto \ket{x + a}$,
    /// where the addition is performed
    /// modulo $2^n$, where $n = \texttt{Length(target!)}$.
    ///
    /// # Input
    /// ## target
    /// A quantum register encoding an unsigned integer using little-endian
    /// encoding in the dual (QFT) basis.
    /// ## increment
    /// The integer by which the `target` is incremented by.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.IncrementByInteger
    ///
    /// # References
    /// - [ *Thomas G. Draper*,
    ///      arXiv:quant-ph/0008033](https://arxiv.org/pdf/quant-ph/0008033v1.pdf)
    ///
    /// # Remarks
    /// Note that we have simplified the circuit because the increment is a
    /// classical constant, not a quantum register.
    ///
    /// See the figure on
    /// [ Page 6 of arXiv:quant-ph/0008033v1 ](https://arxiv.org/pdf/quant-ph/0008033.pdf#page=6)
    /// for the circuit diagram and explanation.
    operation IncrementPhaseByInteger(increment : Int, target : PhaseLittleEndian) : Unit is Adj + Ctl {
        let d = Length(target!);

        for j in 0 .. d - 1 {
            //  Use Microsoft.Quantum.Intrinsic.R1Frac
            R1Frac(increment, (d - 1) - j, (target!)[j]);
        }
    }

    /// # Summary
    /// Increments an unsigned quantum register by a classical integer,
    /// using phase rotations.
    ///
    /// # Description
    /// Suppose that `target` encodes an unsigned integer $x$ in a little-endian
    /// encoding and that `increment` is equal to $a$.
    /// Then, this operation implements the unitary $\ket{x} \mapsto \ket{x + a}$,
    /// where the addition is performed
    /// modulo $2^n$, where $n = \texttt{Length(target!)}$.
    ///
    /// # Input
    /// ## target
    /// A quantum register encoding an unsigned integer using little-endian
    /// encoding.
    /// ## increment
    /// The integer by which the `target` is incremented by.
    operation IncrementByInteger(increment : Int, target : LittleEndian) : Unit is Adj + Ctl {
        ApplyPhaseLEOperationOnLECA(IncrementPhaseByInteger(increment, _), target);
    }

}
