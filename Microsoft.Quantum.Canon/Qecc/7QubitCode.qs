// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Private operation used to implement both the Steane code encoder and decoder.
    ///
    /// # Input
    /// ## data
    /// an array holding 1 qubit which is the input qubit.
    /// ## scratch
    /// an array holding 6 qubits which add redundancy.
    operation SteaneCodeEncoderImpl(data : Qubit[], scratch : Qubit[])  : ()
    {
        body {
            H( scratch[0] );
            H( scratch[2] );
            H( scratch[5] );
            CNOT( data[0],   scratch[4] );
            CNOT( scratch[5], scratch[1] );
            CNOT( scratch[5], scratch[3] );
            CNOT( scratch[1], data[0]   );
            CNOT( scratch[2], scratch[4] );
            CNOT( scratch[0], scratch[4] );
            CNOT( scratch[4], scratch[5] );
            CNOT( scratch[2], scratch[3] );
            CNOT( scratch[0], scratch[1] );
        }
        adjoint auto
    }

    /// # Summary
    /// Decoder for the ⟦7, 1, 3⟧ Steane quantum code.
    /// 
    /// # Remarks
    /// The chosen decoder uses the CSS code property of the ⟦7, 1, 3⟧ Steane code, i.e., it corrects X errors 
    /// and Z errors separately. A property of the code is that the location of the X, respectively, Z correction 
    /// to be applied is the 3-bit encoding of the X, repsectively, Z syndrome when considered an integer. For more 
    /// information, see D. Gottesman, "Stabilizer Codes and Quantum Error Correction," Ph.D. Thesis, Caltech, 1997; 
    /// https://arxiv.org/abs/quant-ph/9705052
    function SteaneCodeRecoveryX( syndrome : Syndrome)  : Pauli[]
    {
        return EmbedPauli(PauliX, ResultAsInt(syndrome), 7);
    }

    function SteaneCodeRecoveryZ( syndrome : Syndrome)  : Pauli[]
    {
        return EmbedPauli(PauliZ, ResultAsInt(syndrome), 7);
    }

    function SteaneCodeRecoveryFns() : (RecoveryFn, RecoveryFn) {
        return (RecoveryFn(SteaneCodeRecoveryX), RecoveryFn(SteaneCodeRecoveryZ));
    }

    /// # Summary
    /// Encodes into the ⟦7, 1, 3⟧ Steane quantum code.
    operation SteaneCodeEncoder(physRegister : Qubit[], auxQubits : Qubit[])  : LogicalRegister
    {
        body {
            SteaneCodeEncoderImpl(physRegister, auxQubits);

            let logicalRegister = LogicalRegister(physRegister + auxQubits);
            return logicalRegister;
        }
    }

    operation SteaneCodeDecoder( logicalRegister : LogicalRegister)  : (Qubit[], Qubit[])
    {
        body {
            let physRegister = [logicalRegister[0]];
            let auxQubits = logicalRegister[1..6];

            (Adjoint SteaneCodeEncoderImpl)(physRegister, auxQubits);

            return (physRegister, auxQubits);
        }
    }

    /// # Summary
    /// Returns a CSS value representing the ⟦7, 1, 3⟧ Steane code encoder and
    /// decoder with in-place syndrome measurement.
    operation  SteaneCode()  : CSS
    {
        body {
            let e = EncodeOp(SteaneCodeEncoder);
            let d = DecodeOp(SteaneCodeDecoder);
            let sx = SyndromeMeasOp(MeasureStabilizerGenerators(
                        [ [ PauliX; PauliI; PauliX; PauliI; PauliX; PauliI; PauliX ];
                        [ PauliI; PauliX; PauliX; PauliI; PauliI; PauliX; PauliX ];
                        [ PauliI; PauliI; PauliI; PauliX; PauliX; PauliX; PauliX ] ],
                        _, MeasureWithScratch)
                    );
            let sz = SyndromeMeasOp(MeasureStabilizerGenerators(
                        [ [ PauliZ; PauliI; PauliZ; PauliI; PauliZ; PauliI; PauliZ ];
                        [ PauliI; PauliZ; PauliZ; PauliI; PauliI; PauliZ; PauliZ ];
                        [ PauliI; PauliI; PauliI; PauliZ; PauliZ; PauliZ; PauliZ ] ],                    
                        _, MeasureWithScratch)
                    );
            let code = CSS(e, d, sx, sz);
            return code;
        }
    }


}
