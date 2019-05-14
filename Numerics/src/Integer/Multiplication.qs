// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
	open Microsoft.Quantum.Arrays;
	open Microsoft.Quantum.Intrinsic;
	open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Multiply integer `xs` by integer `ys` and store the result in `result`,
    /// which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// $n$-bit multiplicand (LittleEndian)
    /// ## ys
    /// $n$-bit multiplier (LittleEndian)
    /// ## result
    /// $2n$-bit result (LittleEndian), must be in state $\ket{0}$ initially.
    ///
    /// # Remarks
    /// Uses a standard shift-and-add approach to implement the multiplication.
    /// The controlled version was improved by copying out $x_i$ to an ancilla
    /// qubit conditioned on the control qubits, and then controlling the
    /// addition on the ancilla qubit.
    operation IntegerMultiplication (xs: LittleEndian, ys: LittleEndian,
                                     result: LittleEndian) : Unit {
        body (...) {
            let n = Length(xs!);

            EqualityFactI(n, Length(ys!), "Integer multiplication requires
                           equally-sized registers xs and ys.");
            EqualityFactI(2 * n, Length(result!), "Integer multiplication
                            requires a 2n-bit result registers.");
            AssertAllZero(result!);

			for (i in 0..n-1) {
                (Controlled IntegerAddition) ([xs![i]], (ys,
                    LittleEndian(result![i..i+n])));
            }
        }
        controlled (controls, ...) {
            let n = Length(xs!);

            EqualityFactI(n, Length(ys!), "Integer multiplication requires
                           equally-sized registers xs and ys.");
            EqualityFactI(2 * n, Length(result!), "Integer multiplication
                            requires a 2n-bit result registers.");
            AssertAllZero(result!);

            using (anc = Qubit()) {
			    for (i in 0..n-1) {
                    (Controlled CNOT) (controls, (xs![i], anc));
                    (Controlled IntegerAddition) ([anc], (ys,
                        LittleEndian(result![i..i+n])));
                    (Controlled CNOT) (controls, (xs![i], anc));
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Computes the square of the integer `xs` into `result`,
    /// which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// $n$-bit number to square (LittleEndian)
    /// ## result
    /// $2n$-bit result (LittleEndian), must be in state $\ket{0}$ initially.
    ///
    /// # Remarks
    /// Uses a standard shift-and-add approach to compute the square. Saves
    /// $n-1$ qubits compared to the straight-forward solution which first
    /// copies out xs before applying a regular multiplier and then undoing
    /// the copy operation.
    operation IntegerSquare (xs: LittleEndian, result: LittleEndian) : Unit {
        body (...) {
            (Controlled IntegerSquare) (new Qubit[0], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!);

            EqualityFactI(2 * n, Length(result!), "Integer multiplication
                            requires a 2n-bit result registers.");
            AssertAllZero(result!);

            using (anc = Qubit()) {
			    for (i in 0..n-1) {
                    (Controlled CNOT) (controls, (xs![i], anc));
                    (Controlled IntegerAddition) ([anc], (xs,
                        LittleEndian(result![i..i+n])));
                    (Controlled CNOT) (controls, (xs![i], anc));
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Multiply signed integer `xs` by signed integer `ys` and store
    /// the result in `result`, which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// n-bit multiplicand (SignedLittleEndian)
    /// ## ys
    /// n-bit multiplier (SignedLittleEndian)
    /// ## result
    /// 2n-bit result (SignedLittleEndian), must be in state $\ket{0}$
    /// initially.
    operation SignedIntegerMultiplication (xs: SignedLittleEndian,
                                           ys: SignedLittleEndian,
                                           result: SignedLittleEndian): Unit {
        body (...) {
            (Controlled SignedIntegerMultiplication) (new Qubit[0],
                                                      (xs, ys, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!!);
            using ((signx, signy) = (Qubit(), Qubit())) {
                CNOT(Tail(xs!!), signx);
                CNOT(Tail(ys!!), signy);
                (Controlled IntegerInversion2s)([signx], xs);
                (Controlled IntegerInversion2s)([signy], ys);

                (Controlled IntegerMultiplication) (controls,
                                                    (xs!, ys!, result!));
                CNOT(signx, signy);
                // No controls required since `result` will still be zero
                // if we did not perform the multiplication above.
                (Controlled IntegerInversion2s)([signy], result);
                CNOT(signx, signy);

                (Controlled Adjoint IntegerInversion2s)([signx], xs);
                (Controlled Adjoint IntegerInversion2s)([signy], ys);
                CNOT(Tail(xs!!), signx);
                CNOT(Tail(ys!!), signy);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Square signed integer `xs` and store
    /// the result in `result`, which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// n-bit integer to square (SignedLittleEndian)
    /// ## result
    /// 2n-bit result (SignedLittleEndian), must be in state $\ket{0}$
    /// initially.
    ///
    /// # Remarks
    /// The implementation relies on IntegerSquare.
    operation SignedIntegerSquare (xs: SignedLittleEndian,
                                   result: SignedLittleEndian): Unit {
        body (...) {
            (Controlled SignedIntegerSquare) (new Qubit[0], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!!);
            using ((signx, signy) = (Qubit(), Qubit())) {
                CNOT(Tail(xs!!), signx);
                (Controlled IntegerInversion2s)([signx], xs);

                (Controlled IntegerSquare) (controls, (xs!, result!));

                (Controlled Adjoint IntegerInversion2s)([signx], xs);
                CNOT(Tail(xs!!), signx);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}