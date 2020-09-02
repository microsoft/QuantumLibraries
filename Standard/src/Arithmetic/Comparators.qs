// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;

    /// # Summary
    /// This operation tests if an integer represented by a register of qubits
    /// is greater than another integer, applying an XOR of the result onto an
    /// output qubit.
    ///
    /// # Description
    /// Given two integers `x` and `y` stored in equal-size qubit registers,
    /// this operation checks if they satisfy `x > y`. If true, 1 is
    /// XORed into an output qubit. Otherwise, 0 is XORed into an output qubit.
    /// In other words, this operation can be represented by the unitary
    /// $$
    /// \begin{align}
    ///     U\ket{x}\ket{y}\ket{z} = \ket{x}\ket{y}\ket{z\oplus (x>y)}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## x
    /// First number to be compared stored in `LittleEndian` format in a qubit register.
    /// ## y
    /// Second number to be compared stored in `LittleEndian` format in a qubit register.
    /// ## output
    /// Qubit that stores the result of the comparison $x>y$.
    ///
    /// # References
    /// - A new quantum ripple-carry addition circuit
    ///   Steven A. Cuccaro, Thomas G. Draper, Samuel A. Kutin, David Petrie Moulton
    ///   https://arxiv.org/abs/quant-ph/0410184
    operation CompareUsingRippleCarry(x : LittleEndian, y : LittleEndian, output : Qubit)
    : Unit is Adj + Ctl {
        if (Length(x!) != Length(y!)) {
            fail "Size of integer registers must be equal.";
        }

        using (auxiliary = Qubit()) {
            within {
                let nQubitsX = Length(x!);

                // Take 2's complement
                ApplyToEachCA(X, x! + [auxiliary]);

                ApplyMajorityInPlace(x![0], [y![0], auxiliary]);
                ApplyToEachCA(MAJ, Zip3(Most(x!), Rest(y!), Rest(x!)));
            } apply {
                X(output);
                CNOT(Tail(x!), output);
            }
        }
    }

    /// # Summary
    /// This operation tests if an integer represented by a register of qubits is
    /// less than a big integer provided as a constant.
    ///
    /// # Description
    /// Given two integers `x` and `c`, `x` stored in a qubit register, and `c` being
    /// a big integer constant, this operation checks if they satisfy `x < c`.  If true,
    /// the output qubit is changed to state $\ket 1$.  The output qubit is assumed to
    /// be in state $\ket 0$ when the operation is being called.
    ///
    /// # Input
    /// ## c
    /// Non-negative constant number to compared to
    /// ## x
    /// Number in qubit register to compare
    /// ## output
    /// Qubit that stores the result of comparison (must be initialized to $\ket 0$)
    ///
    /// # Remark
    /// This operation applies several optimizations in addition to the construction
    /// described in the original work.  Special cases are applied if `c` is $0$,
    /// `c` is $2^n$, or `c` is $2^{n-1}$, where $n$ is the number of bits in `x`.
    /// Qubits and AND gates are saved for each trailing `0` in the bit representation of
    /// `c`.  Further, one AND gate is saved for the least significant bit in `c`, which
    /// is 1, after the trailing `0`s have been removed.  Further, all qubits associated
    /// to constant inputs in the construction of the original work are propagated and
    /// not allocated in this implementation.
    ///
    /// # References
    /// - Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization
    ///   Dominic W. Berry, Craig Gidney, Mario Motta, Jarrod R. McClean, Ryan Babbush
    ///   Quantum 3, 208 (2019)
    ///   https://arxiv.org/abs/1902.02134v4
    operation CompareLessThanConstantUsingRippleCarry(c : BigInt, x : LittleEndian, output : Qubit)
    : Unit is Adj+Ctl {
        body (...) {
            let bitwidth = Length(x!);
            AssertAllZero([output]);
            Fact(c >= 0L, $"Constant input `c` must not be negative, but was {c}");

            if (c == 0L) {
                // do nothing; output stays 0
            } elif (c >= (1L <<< bitwidth)) {
                X(output);
            } elif (c == (1L <<< (bitwidth - 1))) {
                CNOT(Tail(x!), output);
                X(output);
            } else {
                let bits = BigIntAsBoolArray(c);
                let l = IndexOf(Identity<Bool>, bits);
                using (tmpAnd = Qubit[bitwidth - 2 - l]) {
                    let tmpCarry = x![l..l] + tmpAnd;
                    within {
                        ApplyPauliFromBitString(PauliX, true, bits, x!);
                        ApplyToEachA(ApplyAndOrStep, Zip4(bits[l + 1...], Most(tmpCarry), x![l + 1...], Rest(tmpCarry)));
                    } apply {
                        ApplyAndOrStep(bits[bitwidth - 1], Tail(tmpCarry), Tail(x!), output);
                    }
                }
            }
        }
        adjoint self;

        controlled (controls, ...) {
            using (q = Qubit()) {
                within {
                    CompareLessThanConstantUsingRippleCarry(c, x, q);
                } apply {
                    if (Length(controls) == 1) {
                        ApplyAnd(Head(controls), q, output);
                    } else {
                        Controlled X(controls + [q], output);
                    }
                }
            }
        }
        adjoint controlled self;
    }

    /// # Summary
    /// Applies to input-transformed AND or OR gate
    ///
    /// # Description
    /// Given a $\ket 0$ initialized `target`, applies
    /// $a \land \bar b$, if `isOr` is `false`, and
    /// $a \lor b$, if `isOr` is `true`
    ///
    /// # Remark
    /// The OR gate is realized using `ApplyAnd`
    internal operation ApplyAndOrStep(isOr : Bool, a : Qubit, b : Qubit, target : Qubit)
    : Unit is Adj {
        within {
            ApplyIfA(X, isOr, a);
            X(b);
        } apply {
            ApplyAnd(a, b, target);
            ApplyIfA(X, isOr, target);
        }
    }

}
