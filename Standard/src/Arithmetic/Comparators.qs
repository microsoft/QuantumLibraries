// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;

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

    operation LessThanConstantUsingRippleCarry(c : BigInt, x : LittleEndian, output : Qubit)
    : Unit {
        let bitwidth = Length(x!);
        AssertAllZero([output]);

        if (c == 0L) {
            // do nothing; output stays 0
        } elif (c >= (1L <<< bitwidth)) {
            X(output);
        } elif (c == (1L <<< (bitwidth - 1))) {
            CNOT(Tail(x!), output);
            X(output);
        } else {
            let l = TrailingZeroes(c);
            using ((tmpConstants, tmpAnd) = (Qubit[bitwidth - 1 - l], Qubit[bitwidth - 2 - l])) {
                let tmpCarry = x![l..l] + tmpAnd;
                within {
                    ApplyXorInPlaceL(c >>> l + 1, LittleEndian(tmpConstants));
                    ApplyXorInPlaceL(c, x);
                    for (i in 0..bitwidth - l - 2) {
                        if (i > 0) {
                            ApplyAnd(tmpConstants[i - 1], x![i + l], tmpCarry[i]);
                            CNOT(tmpCarry[i - 1], tmpCarry[i]);
                        }
                        CNOT(tmpCarry[i], tmpConstants[i]);
                    }
                } apply {
                    ApplyAnd(Tail(tmpConstants), Tail(x!), output);
                    CNOT(Tail(tmpCarry), output);
                }
            }
        }
    }

    internal function TrailingZeroes(number : BigInt) : Int {
        mutable zeroes = 0;
        mutable copy = number;
        while (copy % 2L == 0L) {
            set zeroes += 1;
            set copy /= 2L;
        }
        return zeroes;
    }

}
