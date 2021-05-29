// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Inverts a given integer modulo 2's complement.
    ///
    /// # Input
    /// ## xs
    /// n-bit signed integer (SignedLittleEndian), will be inverted modulo
    /// 2's complement.
    operation Invert2sSI (xs: SignedLittleEndian) : Unit is Adj + Ctl {
        body (...) {
            Controlled Invert2sSI([], xs);
        }
        controlled (controls, ...) {
            ApplyToEachCA((Controlled X)(controls, _), xs!!);

            use aux = Qubit[Length(xs!!)];
            within {
                Controlled X(controls, aux[0]);
            } apply {
                AddI(LittleEndian(aux), xs!);
            }
        }
    }

    /// # Summary
    /// Computes the reciprocal 1/x for an unsigned integer x
    /// using integer division. The result, interpreted as an integer,
    /// will be `floor(2^(2*n-1) / x)`.
    ///
    /// # Input
    /// ## xs
    /// n-bit unsigned integer
    /// ## result
    /// 2n-bit output, must be in $\ket{0}$ initially.
    ///
    /// # Remarks
    /// For the input x=0, the output will be all-ones.
    operation ComputeReciprocalI (xs: LittleEndian,
                                  result: LittleEndian) : Unit is Adj + Ctl {
        body (...) {
            Controlled ComputeReciprocalI([], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!);
            EqualityFactI(Length(result!), 2*n,
                          "Result register must contain 2n qubits.");
            AssertAllZero(result!);
            use lhs = Qubit[2 * n];
            use padding = Qubit[n];
            let paddedxs = LittleEndian(xs! + padding);
            X(Tail(lhs)); // initialize left-hand side to 2^{2n-1}
            // ... and divide:
            (Controlled DivideI) (controls,
                (LittleEndian(lhs), paddedxs, result));
            // uncompute lhs
            for i in 0..2 * n - 1 {
                (Controlled AddI) ([result![i]],
                    (LittleEndian(paddedxs![0..2*n-1-i]),
                        LittleEndian(lhs[i..2*n-1])));
            }
            X(Tail(lhs));
        }
    }
}
