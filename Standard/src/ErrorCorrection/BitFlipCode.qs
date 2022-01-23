// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Private operation used to implement both the bit flip encoder and decoder.
    ///
    /// # Remarks
    /// Note that this encoder can make use of in-place coherent recovery,
    /// in which case it will "cause" the error described
    /// by the initial state of `auxQubits`.
    /// In particular, if `auxQubits` are initially in the state $\ket{10}$, this
    /// will cause an $X_1$ error on the encoded qubit.
    ///
    /// # References
    /// - doi:10.1103/PhysRevA.85.044302
    internal operation ApplyBitFlipEncoder(coherentRecovery : Bool, data : Qubit[], scratch : Qubit[])
    : Unit is Adj {
        if coherentRecovery {
            Controlled X(scratch, data[0]);
        }

        Controlled X(data, scratch[0]);
        Controlled X(data, scratch[1]);
    }


    /// # Summary
    /// Encodes into the [3, 1, 3] / ⟦3, 1, 1⟧ bit-flip code.
    ///
    /// # Input
    /// ## physRegister
    /// A register of physical qubits representing the data to be protected.
    /// ## auxQubits
    /// A register of auxiliary qubits initially in the $\ket{00}$ state to be
    /// used in encoding the data to be protected.
    ///
    /// # Output
    /// The physical and auxiliary qubits used in encoding, represented as a
    /// logical register.
    ///
    /// # See Also
    /// - Microsoft.Quantum.ErrorCorrection.LogicalRegister
    operation EncodeIntoBitFlipCode(physRegister : Qubit[], auxQubits : Qubit[])
    : LogicalRegister {
        ApplyBitFlipEncoder(false, physRegister, auxQubits);
        let logicalRegister = LogicalRegister(physRegister + auxQubits);
        return logicalRegister;
    }


    /// # Summary
    /// Decodes from the [3, 1, 3] / ⟦3, 1, 1⟧ bit-flip code.
    ///
    /// # Input
    /// ## logicalRegister
    /// A code block of the bit-flip code.
    ///
    /// # Output
    /// A tuple of the data encoded into the logical register, and the auxiliary
    /// qubits used to represent the syndrome.
    ///
    /// # See Also
    /// - Microsoft.Quantum.ErrorCorrection.LogicalRegister
    /// - Microsoft.Quantum.ErrorCorrection.EncodeIntoBitFlipCode
    operation DecodeFromBitFlipCode(logicalRegister : LogicalRegister)
    : (Qubit[], Qubit[]) {
        let physRegister = [(logicalRegister!)[0]];
        let auxQubits = (logicalRegister!)[1 .. 2];
        Adjoint ApplyBitFlipEncoder(false, physRegister, auxQubits);
        return (physRegister, auxQubits);
    }

    /// # Summary
    /// Returns a QECC value representing the ⟦3, 1, 1⟧ bit flip code encoder and
    /// decoder with in-place syndrome measurement.
    ///
    /// # Output
    /// Returns an implementation of a quantum error correction code by
    /// specifying a `QECC` type.
    function BitFlipCode () : QECC {
        let e = EncodeOp(EncodeIntoBitFlipCode);
        let d = DecodeOp(DecodeFromBitFlipCode);
        let s = SyndromeMeasOp(MeasureStabilizerGenerators([[PauliZ, PauliZ, PauliI], [PauliI, PauliZ, PauliZ]], _, MeasureWithScratch));
        let code = QECC(e, d, s);
        return code;
    }

    /// # Summary
    /// Function for recovery Pauli operations for given syndrome measurement
    /// by table lookup for the ⟦3, 1, 1⟧ bit flip code.
    ///
    /// # Output
    /// Function of type `RecoveryFn` that takes a syndrome measurement
    /// `Result[]` and returns the `Pauli[]` operations that corrects the
    /// detected error.
    ///
    /// # See Also
    /// - Microsoft.Quantum.ErrorCorrection.RecoveryFn
    function BitFlipRecoveryFn () : RecoveryFn {
        return TableLookupRecovery([
            [PauliI, PauliI, PauliI],
            [PauliX, PauliI, PauliI],
            [PauliI, PauliI, PauliX],
            [PauliI, PauliX, PauliI]
        ]);
    }

}
