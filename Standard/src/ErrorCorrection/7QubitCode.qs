// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Private operation used to implement both the Steane code encoder and decoder.
    ///
    /// # Input
    /// ## data
    /// an array holding 1 qubit which is the input qubit.
    /// ## scratch
    /// an array holding 6 qubits which add redundancy.
    operation SteaneCodeEncoderImpl (data : Qubit[], scratch : Qubit[]) : Unit {
        body (...) {
            H(scratch[0]);
            H(scratch[2]);
            H(scratch[5]);
            CNOT(data[0], scratch[4]);
            CNOT(scratch[5], scratch[1]);
            CNOT(scratch[5], scratch[3]);
            CNOT(scratch[1], data[0]);
            CNOT(scratch[2], scratch[4]);
            CNOT(scratch[0], scratch[4]);
            CNOT(scratch[4], scratch[5]);
            CNOT(scratch[2], scratch[3]);
            CNOT(scratch[0], scratch[1]);
        }

        adjoint invert;
    }

    /// # Summary
    /// Decoder for the X-part of the stabilizer group of the ⟦7, 1, 3⟧ Steane quantum code.
    ///
    /// # Input
    /// ## syndrome
    /// A syndrome array obtained from measuring the X-part of the stabilizer.
    ///
    /// # Output
    /// An array of Pauli operations which, when applied to the encoded quantum system
    /// corrects the error corresponding to `syndrome`.
    ///
    /// # Remarks
    /// The chosen decoder uses the CSS code property of the ⟦7, 1, 3⟧ Steane code, i.e., it corrects X errors
    /// and Z errors separately. A property of the code is that the location of the X, respectively, Z correction
    /// to be applied is the 3-bit encoding of the X, respectively, Z syndrome when considered an integer.
    ///
    /// # See Also
    /// - SteaneCodeRecoveryX
    /// - SteaneCodeRecoveryFns
    ///
    /// # References
    /// - D. Gottesman, "Stabilizer Codes and Quantum Error Correction," Ph.D. Thesis, Caltech, 1997;
    /// https://arxiv.org/abs/quant-ph/9705052
    function SteaneCodeRecoveryX (syndrome : Syndrome) : Pauli[] {
        let idxQubit = ResultArrayAsInt(syndrome!);
        return idxQubit == 0
               ? ConstantArray(7, PauliI)
               | EmbedPauli(PauliZ, idxQubit - 1, 7);
    }

    /// # Summary
    /// Decoder for the Z-part of the stabilizer group of the ⟦7, 1, 3⟧ Steane quantum code.
    ///
    /// # See Also
    /// - SteaneCodeRecoveryX
    /// - SteaneCodeRecoveryFns
    function SteaneCodeRecoveryZ (syndrome : Syndrome) : Pauli[] {
        let idxQubit = ResultArrayAsInt(syndrome!);
        return idxQubit == 0
               ? ConstantArray(7, PauliI)
               | EmbedPauli(PauliX, idxQubit - 1, 7);
    }

    /// # Summary
    /// Decoder for combined X- and Z-parts of the stabilizer group of the
    /// ⟦7, 1, 3⟧ Steane quantum code.
    ///
    /// # Output
    /// Tuple of functions of type `RecoveryFn` that takes a syndrome
    /// measurement `Result[]` and returns the `Pauli[]` operations that
    /// corrects the detected error.
    ///
    /// # See Also
    /// - SteaneCodeRecoveryX
    /// - SteaneCodeRecoveryZ
    function SteaneCodeRecoveryFns () : (RecoveryFn, RecoveryFn)
    {
        return (RecoveryFn(SteaneCodeRecoveryX), RecoveryFn(SteaneCodeRecoveryZ));
    }
    
    
    /// # Summary
    /// An encoding operation that maps an unencoded quantum register to an encoded quantum register
    /// under the ⟦7, 1, 3⟧ Steane quantum code.
    ///
    /// # Input
    /// ## physRegister
    /// A qubit register which holds the an unencoded quantum state
    /// ## auxQubits
    /// A qubit register which is initially zero and which gets added to the quantum
    /// system so that an encoding operation can be performed
    ///
    /// # Output
    /// A quantum register holding the state after the Steane encoder has been applied
    ///
    /// # See Also
    /// - LogicalRegister
    /// - SteaneCodeDecoder
    operation EncodeIntoSteaneCode(physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister
    {
        SteaneCodeEncoderImpl(physRegister, auxQubits);
        let logicalRegister = LogicalRegister(physRegister + auxQubits);
        return logicalRegister;
    }
    
    
    /// # Summary
    /// An inverse encoding operation that maps an unencoded quantum register to an encoded quantum
    /// register under the ⟦7, 1, 3⟧ Steane quantum code.
    ///
    /// # Input
    /// ## logicalRegister
    /// An array of qubits representing the encoded 5-qubit code logical state.
    ///
    /// # Output
    /// A qubit array of length 1 representing the unencoded state in the
    /// first parameter, together with auxiliary qubits in the second parameter.
    ///
    /// # Remarks
    /// The chosen decoder uses the CSS code property of the ⟦7, 1, 3⟧ Steane code, i.e., it corrects X errors
    /// and Z errors separately. A property of the code is that the location of the X, respectively, Z correction
    /// to be applied is the 3-bit encoding of the X, respectively, Z syndrome when considered an integer.
    ///
    /// # See Also
    /// - SteaneCodeEncoder
    /// - SteaneCodeDecoder
    /// - LogicalRegister
    ///
    /// # References
    /// - D. Gottesman, "Stabilizer Codes and Quantum Error Correction," Ph.D. Thesis, Caltech, 1997;
    ///   https://arxiv.org/abs/quant-ph/9705052
    operation DecodeFromSteaneCode(logicalRegister : LogicalRegister) : (Qubit[], Qubit[])
    {
        let physRegister = [(logicalRegister!)[0]];
        let auxQubits = (logicalRegister!)[1 .. 6];
        Adjoint SteaneCodeEncoderImpl(physRegister, auxQubits);
        return (physRegister, auxQubits);
    }

    /// # Summary
    /// Returns a CSS value representing the ⟦7, 1, 3⟧ Steane code encoder and
    /// decoder with in-place syndrome measurement.
    ///
    /// # Output
    /// An object of CSS type which collects all relevant data to perform encoding and
    /// error correction for the ⟦7, 1, 3⟧ Steane code.
    ///
    /// # Remarks
    /// This code was found in the following paper:
    /// - A. Steane, "Multiple Particle Interference and Quantum Error Correction", Proc. Roy. Soc. Lond. A452 (1996) pp. 2551; https://arxiv.org/abs/quant-ph/9601029.
    function SteaneCode () : CSS {
        let e = EncodeOp(EncodeIntoSteaneCode);
        let d = DecodeOp(DecodeFromSteaneCode);
        let xg = 
            [ 
                [PauliX, PauliI, PauliX, PauliI, PauliX, PauliI, PauliX],
                [PauliI, PauliX, PauliX, PauliI, PauliI, PauliX, PauliX],
                [PauliI, PauliI, PauliI, PauliX, PauliX, PauliX, PauliX]
            ];
        let zg = 
            [ 
                [PauliZ, PauliI, PauliZ, PauliI, PauliZ, PauliI, PauliZ],
                [PauliI, PauliZ, PauliZ, PauliI, PauliI, PauliZ, PauliZ],
                [PauliI, PauliI, PauliI, PauliZ, PauliZ, PauliZ, PauliZ] 
            ];
        let sx = SyndromeMeasOp(MeasureStabilizerGenerators(xg, _, MeasureWithScratch));
        let sz = SyndromeMeasOp(MeasureStabilizerGenerators(zg, _, MeasureWithScratch));
        let code = CSS(e, d, sx, sz);
        return code;
    }

}


