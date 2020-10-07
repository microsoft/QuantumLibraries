// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies a bitwise-XOR operation between a classical integer and an
    /// integer represented by a register of qubits.
    ///
    /// # Description
    /// Applies `X` operations to qubits in a little-endian register based on
    /// 1 bits in an integer.
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
    operation ApplyXorInPlace(value : Int, target : LittleEndian)
    : Unit is Adj + Ctl {
        ApplyToEachCA(
            CControlledCA(X),
            Zipped(IntAsBoolArray(value, Length(target!)), target!)
        );
    }

    /// # Summary
    /// Applies the three-qubit majority operation in-place on a register of
    /// qubits.
    ///
    /// # Description
    /// This operation computes the majority function in-place on 3 qubits.
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
    operation ApplyMajorityInPlace(output: Qubit, input: Qubit[])
    : Unit is Adj + Ctl {
        if (Length(input) == 2) {
            MAJ(input[0], input[1], output);
        } else {
            fail $"The in-place majority operation on {Length(input)} is qubits not yet implemented.";
        }
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
    operation CopyMostSignificantBit (from : LittleEndian, target : Qubit) : Unit {
        body (...) {
            let mostSingificantQubit = Tail(from!);
            CNOT(mostSingificantQubit, target);
        }

        adjoint invert;
    }

}


