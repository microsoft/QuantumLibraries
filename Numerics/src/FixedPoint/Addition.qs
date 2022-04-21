// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Adds a classical constant to a quantum fixed-point number.
    ///
    /// # Input
    /// ## constant
    /// Constant to add to the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number to which the constant will
    /// be added.
    operation AddConstantFxP(constant : Double, fp : FixedPoint) : Unit is Adj + Ctl {
        let (px, xs) = fp!;
        let n = Length(xs);
        use ys = Qubit[n];
        let tmpFp = FixedPoint(px, ys);
        ApplyWithCA(PrepareFxP(constant, _), AddFxP(_, fp), tmpFp);
    }

    /// # Summary
    /// Adds two fixed-point numbers stored in quantum registers.
    ///
    /// # Description
    /// Given two fixed-point registers respectively in states $\ket{f_1}$ and $\ket{f_2}$,
    /// performs the operation $\ket{f_1} \ket{f_2} \mapsto \ket{f_1} \ket{f_1 + f_2}$.
    ///
    /// # Input
    /// ## fp1
    /// First fixed-point number
    /// ## fp2
    /// Second fixed-point number, will be updated to contain the sum of the
    /// two inputs.
    ///
    /// # Remarks
    /// The current implementation requires the two fixed-point numbers
    /// to have the same point position counting from the least-significant
    /// bit, i.e., $n_i$ and $p_i$ must be equal.
    operation AddFxP(fp1 : FixedPoint, fp2 : FixedPoint) : Unit is Adj + Ctl {
        let (px, xs) = fp1!;
        let (py, ys) = fp2!;

        IdenticalPointPosFactFxP([fp1, fp2]);

        AddI(LittleEndian(xs), LittleEndian(ys));
    }

    /// # Summary
    /// Computes the additive inverse of `fp`.
    ///
    /// # Input
    /// ## fp
    /// Fixed-point number to invert.
    ///
    /// # Remarks
    /// Numerical inaccuracies may occur depending on the
    /// bit-precision of the fixed-point number.
    ///
    /// # See also
    /// - Microsoft.Quantum.Arithmetic.SubtractFxP
    operation InvertFxP(fp: FixedPoint) : Unit is Adj + Ctl {
        let (_, reg) = fp!;
        Invert2sSI(SignedLittleEndian(LittleEndian(reg)));
    }

    /// # Summary
    /// Computes `minuend - subtrahend` and stores the difference in `minuend`.
    ///
    /// # Input
    /// ## subtrahend
    /// The subtrahend of the subtraction - the number to be subtracted.
    /// ## minuend
    /// The minuend of the subtraction - the number from which the other is subtracted.
    ///
    /// # Remarks
    /// Computes the difference by inverting `subtrahend` before and after adding
    /// it to `minuend`.
    ///
    /// # See also
    /// - Microsoft.Quantum.Arithmetic.AddFxP
    /// - Microsoft.Quantum.Arithmetic.InvertFxP
    operation SubtractFxP(subtrahend : FixedPoint, minuend : FixedPoint) : Unit is Adj + Ctl {
        within {
            InvertFxP(subtrahend);
        } apply {
            AddFxP(subtrahend, minuend);
        }
    }
}
