// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Wrapper for addition: Automatically chooses between addition with
    /// carry and without, depending on the register size of `ys`.
    ///
    /// # Input
    /// ## xs
    /// $n$-bit addend (LittleEndian)
    /// ## ys
    /// Addend with at least $n$ qubits (LittleEndian). Will hold the result.
    operation IntegerAddition (xs: LittleEndian, ys: LittleEndian) : Unit {
        body (...) {
            if (Length(xs!) == Length(ys!)) {
                RippleCarryAdderNoCarryTTK(xs, ys);
            }
            elif (Length(ys!) > Length(xs!)) {
                using (qs = Qubit[Length(ys!) - Length(xs!) - 1]){
                    RippleCarryAdderTTK(LittleEndian(xs! + qs),
                                        LittleEndian(Most(ys!)), Tail(ys!));
                }
            }
            else{
                fail "xs must not contain more qubits than ys!";
            }
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Wrapper for integer comparison: `result = x > y`.
    ///
    /// # Input
    /// ## xs
    /// First $n$-bit number
    /// ## ys
    /// Second $n$-bit number
    /// ## result
    /// Will be flipped if $x > y$
    operation IntegerGreaterThan (xs: LittleEndian, ys: LittleEndian,
                                  result: Qubit) : Unit {
        body (...) {
            GreaterThan(xs, ys, result);
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
    }

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
    /// Integer division of integer `xs` by integer `ys`. `xs` will hold the
    /// remainder `xs - floor(xs/ys) * ys` and `result` will hold
    /// `floor(xs/ys)`.
    ///
    /// # Input
    /// ## xs
    /// n-bit dividend (LittleEndian), will be replaced by the remainder.
    /// ## ys
    /// n-bit divisor (LittleEndian)
    /// ## result
    /// n-bit result (LittleEndian), must be in state $\ket{0}$ initially
    /// and will be replaced by the result of the integer division.
    ///
    /// # Remarks
    /// Uses a standard shift-and-subtract approach to implement the division.
    /// The controlled version is specialized such the subtraction does not
    /// require additional controls.
    operation IntegerDivision (xs: LittleEndian, ys: LittleEndian,
                               result: LittleEndian) : Unit {
        body (...) {
            (Controlled IntegerDivision) (new Qubit[0], (xs, ys, result));
        }
        controlled (controls, ...) {
            let n = Length(result!);

            EqualityFactI(n, Length(ys!), "Integer division requires
                           equally-sized registers ys and result.");
            EqualityFactI(n, Length(xs!), "Integer division
                            requires an n-bit dividend registers.");
            AssertAllZero(result!);

            let xpadded = LittleEndian(xs! + result!);

			for (i in (n-1)..(-1)..0) {
                let xtrunc = LittleEndian(xpadded![i..i+n-1]);
                (Controlled IntegerGreaterThan) (controls,
                                                 (ys, xtrunc, result![i]));
                // if ys > xtrunc, we don't subtract:
                (Controlled X) (controls, result![i]);
                (Controlled Adjoint IntegerAddition) ([result![i]],
                                                      (ys, xtrunc));
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Type of a signed integer stored in little endian (see LittleEndian).
    newtype SignedLittleEndian = LittleEndian;

    /// # Summary
    /// Wrapper for signed integer comparison: `result = xs > ys`.
    ///
    /// # Input
    /// ## xs
    /// First $n$-bit number
    /// ## ys
    /// Second $n$-bit number
    /// ## result
    /// Will be flipped if $xs > ys$
    operation SignedIntegerGreaterThan (xs: SignedLittleEndian,
                                        ys: SignedLittleEndian,
                                        result: Qubit) : Unit {
        body (...) {
            using (tmp = Qubit()) {
                CNOT(Tail(xs!!), tmp);
                CNOT(Tail(ys!!), tmp);
                X(tmp);
                (Controlled IntegerGreaterThan)([tmp], (xs!, ys!, result));
                X(tmp);
                CCNOT(tmp, Tail(ys!!), result);
                CNOT(Tail(xs!!), tmp);
                CNOT(Tail(ys!!), tmp);
            }
        }
        controlled auto;
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

    /// # Summary
    /// Inverts a given integer modulo 2's complement.
    ///
    /// # Input
    /// ## xs
    /// n-bit signed integer (SignedLittleEndian), will be inverted modulo
    /// 2's complement.
    operation IntegerInversion2s (xs: SignedLittleEndian) : Unit {
        body (...) {
            (Controlled IntegerInversion2s) (new Qubit[0], xs);
        }
        controlled (controls, ...) {
            ApplyToEachCA((Controlled X)(controls, _), xs!!);

            using (ancillas = Qubit[Length(xs!!)]) {
                (Controlled X)(controls, ancillas[0]);
                IntegerAddition(LittleEndian(ancillas), xs!);
                (Controlled X)(controls, ancillas[0]);
            }
        }
        adjoint auto;
        adjoint controlled auto;
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
    operation IntegerReciprocal (xs: LittleEndian,
                                 result: LittleEndian) : Unit {
        body (...) {
            (Controlled IntegerReciprocal) (new Qubit[0], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!);
            AssertIntEqual(Length(result!), 2*n,
                           "Result register must contain 2n qubits.");
            AssertAllZero(result!);
            using ((lhs, padding) = (Qubit[2*n], Qubit[n])) {
                let paddedxs = LittleEndian(xs! + padding);
                X(Tail(lhs)); // initialize left-hand side to 2^{2n-1}
                // ... and divide:
                (Controlled IntegerDivision) (controls,
                    (LittleEndian(lhs), paddedxs, result));
                // uncompute lhs
                for (i in 0..2*n-1) {
                    (Controlled IntegerAddition) ([result![i]],
                        (LittleEndian(paddedxs![0..2*n-1-i]),
                         LittleEndian(lhs[i..2*n-1])));
                }
                X(Tail(lhs));
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}