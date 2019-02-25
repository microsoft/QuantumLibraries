// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Applies `X` operations to qubits in a little-endian register based on 1 bits in an integer.
    ///
    /// Let us denote `value` by a and let y be an unsigned integer encoded in `target`,
    /// then `InPlaceXorLE` performs an operation given by the following map:
    /// $\ket{y}\rightarrow \ket{y\oplus a}$ , where $\oplus$ is the bitwise exclusive OR operator.
    ///
    /// # Input
    /// ## value
    /// An integer which is assumed to be non-negative.
    /// ## target
    /// A quantum register which is used to store `value` in little-endian encoding.
    ///
    /// # See Also
    /// - InPlaceXorBE
    operation InPlaceXorLE (value : Int, target : LittleEndian) : Unit {
        body (...) {
            ApplyToEachCA(
                CControlledCA(X),
                Zip(IntAsBoolArray(value, Length(target!)), target!)
            );
        }

        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Summary
    /// Applies `X` operations to qubits in a big-endian register based on 1 bits in an integer.
    ///
    /// Let us denote `value` by a and let y be an unsigned integer encoded in `target`,
    /// then `InPlaceXorBE` performs an operation given by the following map:
    /// $\ket{y}\rightarrow \ket{y\oplus a}$ , where $\oplus$ is the bitwise exclusive OR operator.
    ///
    /// # Input
    /// ## value
    /// An integer which is assumed to be non-negative.
    /// ## target
    /// A quantum register which is used to store `value` in big-endian encoding.
    ///
    /// # See Also
    /// - InPlaceXorLE
    operation InPlaceXorBE(value : Int, target : BigEndian) : Unit {
        body (...) {
            ApplyReversedOpLECA(InPlaceXorLE(value, _), target);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// This computes the Majority function in-place on 3 qubits.
    ///
    /// If we denote output qubit as $z$ and input qubits as $x$ and $y$,
    /// the operation performs the following transformation:
    /// $\ket{xyz} \rightarrow \ket{x \oplus z} \ket{y \oplus z} \ket{\operatorname{MAJ} (x, y, z)}$.
    ///
    /// # Input
    /// ## output
    /// First input qubit. Note that the output is computed in-place
    /// and stored in this qubit.
    /// ## input
    /// Second and third input qubits.
    operation InPlaceMajority(output: Qubit, input: Qubit[]) : Unit {
        body (...) {
            if (Length(input) == 2) {
                MAJ(input[0], input[1], output);
            } else {
                fail $"The in-place majority operation on {Length(input)} is qubits not yet implemented.";
            }
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// This unitary tests if two integers `x` and `y` stored in equal-size qubit registers 
    /// satisfy `x > y`. If true, 1 is XORed into an output
    /// qubit. Otherwise, 0 is XORed into an output qubit.
    ///
    /// In other words, this unitary $U$  satisfies:
    /// $$
    /// \begin{align}
    /// U\ket{x}\ket{y}\ket{z}=\ket{x}\ket{y}\ket{z\oplus (x>y)}.
    /// \end{align}
    /// $$.
    ///
    /// # Input
    /// ## x
    /// First number to be compared stored in `BigEndian` format in a qubit register.
    /// ## y
    /// Second number to be compared stored in `BigEndian` format in a qubit register.
    /// ## output
    /// Qubit that stores the result of the comparison $x>y$.
    ///
    /// # References
    /// - A new quantum ripple-carry addition circuit
    ///   Steven A. Cuccaro, Thomas G. Draper, Samuel A. Kutin, David Petrie Moulton
    ///   https://arxiv.org/abs/quant-ph/0410184
    operation ApplyRippleCarryComparatorBE(x : BigEndian, y : BigEndian, output : Qubit) : Unit {
        body (...) {
            ApplyRippleCarryComparatorLE(BigEndianAsLittleEndian(x), BigEndianAsLittleEndian(y), output);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Copies the most significant bit of a qubit register
    /// `from` representing an unsigned integer into the qubit `target`.
    ///
    /// # Input
    /// ## from
    /// The unsigned integer from which the highest bit is copied from.
    /// the integer is encoded in little-endian format.
    /// ## target
    /// The qubit in which the highest bit is being copied. The bit encoding is
    /// in computational basis.
    ///
    /// # See Also
    /// - LittleEndian
    operation CopyMostSignificantBitLE (from : LittleEndian, target : Qubit) : Unit
    {
        body (...)
        {
            let mostSingificantQubit = Tail(from!);
            CNOT(mostSingificantQubit, target);
        }

        adjoint invert;
    }

}


