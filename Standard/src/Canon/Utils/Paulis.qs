// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Given a single-qubit Pauli operator, applies the corresponding operation
    /// to a single qubit.
    ///
    /// # Input
    /// ## pauli
    /// The Pauli operator to be applied.
    /// ## target
    /// The qubit to which `pauli` is to be applied as an operation.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// ApplyP(PauliX, q);
    /// ```
    /// and
    /// ```qsharp
    /// X(q);
    /// ```
    operation ApplyP(pauli : Pauli, target : Qubit) : Unit is Adj + Ctl {
        if   pauli == PauliX { X(target); }
        elif pauli == PauliY { Y(target); }
        elif pauli == PauliZ { Z(target); }
    }


    /// # Summary
    /// Given a multi-qubit Pauli operator, applies the corresponding operation to
    /// a register.
    ///
    /// # Input
    /// ## pauli
    /// A multi-qubit Pauli operator represented as an array of single-qubit Pauli operators.
    /// ## target
    /// Register to apply the given Pauli operation on.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// ApplyPauli([PauliY, PauliZ, PauliX], target);
    /// ```
    /// and
    /// ```qsharp
    /// Y(target[0]);
    /// Z(target[1]);
    /// X(target[2]);
    /// ```
    operation ApplyPauli(pauli : Pauli[], target : Qubit[]) : Unit is Adj + Ctl {
        ApplyToEachCA(ApplyP, Zipped(pauli, target));
    }


    /// # Summary
    /// Applies a Pauli operator on each qubit in an array if the corresponding
    /// bit of a Boolean array matches a given input.
    ///
    /// # Input
    /// ## pauli
    /// Pauli operator to apply to `qubits[idx]` where `bitsApply == bits[idx]`
    /// ## bitApply
    /// apply Pauli if bit is this value
    /// ## bits
    /// Boolean register specifying which corresponding qubit in `qubits` should be operated on
    /// ## qubits
    /// Quantum register on which to selectively apply the specified Pauli operator
    ///
    /// # Remarks
    /// The Boolean array and the quantum register must be of equal length.
    ///
    /// # Example
    /// The following applies an X operation on qubits 0 and 2, and a Z operation on qubits 1 and 3.
    /// ```qsharp
    /// use qubits = Qubit[4];
    /// let bits = [true, false, true, false];
    /// // Apply when index in `bits` is `true`.
    /// ApplyPauliFromBitString(PauliX, true, bits, qubits);
    /// // Apply when index in `bits` is `false`.
    /// ApplyPauliFromBitString(PauliZ, false, bits, qubits);
    /// ```
    operation ApplyPauliFromBitString(pauli : Pauli, bitApply : Bool, bits : Bool[], qubits : Qubit[])
    : Unit is Adj + Ctl {
        let nBits = Length(bits);

        //FailOn (nbits != Length(qubits), "Number of control bits must be equal to number of control qubits")
        for (bit, qubit) in Zipped(bits, qubits) {
            if bit == bitApply {
                ApplyP(pauli, qubit);
            }
        }
    }


    /// # Summary
    /// Given a single-qubit Pauli operator and the index of a qubit,
    /// returns a multi-qubit Pauli operator with the given single-qubit
    /// operator at that index and `PauliI` at every other index.
    ///
    /// # Input
    /// ## pauli
    /// A single-qubit Pauli operator to be placed at the given location.
    /// ## location
    /// An index such that `output[location] == pauli`, where `output` is
    /// the output of this function.
    /// ## n
    /// Length of the array to be returned.
    ///
    /// # Example
    /// To obtain the array `[PauliI, PauliI, PauliX, PauliI]`:
    /// ```qsharp
    /// EmbedPauli(PauliX, 2, 3);
    /// ```
    function EmbedPauli (pauli : Pauli, location : Int, n : Int) : Pauli[] {
        return ConstantArray(n, PauliI) w/ location <- pauli;
    }

    // FIXME: Remove in favor of something that computes arbitrary
    //        weight Paulis.

    /// # Summary
    /// Returns an array of all weight-1 Pauli operators
    /// on a given number of qubits.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits on which the returned Pauli operators
    /// are defined.
    ///
    /// # Output
    /// An array of multi-qubit Pauli operators, each of which is
    /// represented as an array with length `nQubits`.
    function WeightOnePaulis (nQubits : Int) : Pauli[][] {
        mutable paulis = [[], size = 3 * nQubits];
        let pauliGroup = [PauliX, PauliY, PauliZ];

        for idxQubit in 0 .. nQubits - 1 {
            for idxPauli in IndexRange(pauliGroup) {
                set paulis w/= idxQubit * Length(pauliGroup) + idxPauli <- EmbedPauli(pauliGroup[idxPauli], idxQubit, nQubits);
            }
        }

        return paulis;
    }

}
