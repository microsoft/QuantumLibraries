// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;

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
                ApplyToEachCA(MAJ, Zipped3(Most(x!), Rest(y!), Rest(x!)));
            } apply {
                X(output);
                CNOT(Tail(x!), output);
            }
        }
    }

}
