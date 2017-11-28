// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Private operation used to implement both the 5 qubit encoder and decoder.
    ///
    /// # Input
    /// ## data
    /// an array holding 1 qubit which is the input qubit.
    /// ## scratch
    /// an array holding 4 qubits which add redundancy.
    ///
    /// # Remarks
    /// The particular encoder chosen was taken from the paper V. Kliunchnikov and D. Maslov, "Optimization of Clifford Circuits," 
    /// Phys. Rev. Phys. Rev. A 88, 052307 (2013); https://arxiv.org/abs/1305.0810, Figure 4b) and requires a total of 11 gates.
    operation FiveQubitCodeEncoderImpl(data : Qubit[], scratch : Qubit[])  : ()
    {
        body {
            (Controlled(X))(data, scratch[1]);
            H(data[0]);
            H(scratch[0]);
            (Controlled(X))(data, scratch[2]);
            (Controlled(X))([scratch[0]], data[0]);
            (Controlled(X))(data, scratch[1]);
            (Controlled(X))([scratch[0]], scratch[3]);
            H(scratch[0]);
            H(data[0]);
            (Controlled(X))([scratch[0]], scratch[2]);
            (Controlled(X))(data, scratch[3]);
            // The last X below is to correct the signs of stabilizers.
            // The 5-qubit code is non-CSS, so even if the circuit implements
            // the correct symplectic matrix,
            // it may differ from the desired one by a Pauli correction.
            X(scratch[2]);
        }

        adjoint auto
    }

    /// # Summary
    /// Table lookup decoder for the ⟦5, 1, 3⟧ quantum code.
    ///
    /// # Remarks
    /// By iterating over all errors of weight 1, we obtain a total of 3*5=15 possible non-trivial syndromes. 
    /// Together with the identity, a table of error and corresponding syndrom is built up. For the 5qubit code  
    /// this table is given by: X_1 (0,0,0,1); X_2 (1,0,0,0); X_3 (1,1,0,0); X_4 (0,1,1,0); X_5 (0,0,1,1), 
    /// Z_1 (1,0,1,0); Z_2 (0,1,0,1); Z_3 (0,0,1,0); Z_4 (1,0,0,1); Z_5 (0,1,0,0) with Yi = Xi + Yi. Note that the 
    /// ordering in the table lookup recovery is given by converting the bitvectors to integers (using little endian).
    function  FiveQubitCodeRecoveryFn()  : RecoveryFn
    {
        return TableLookupRecovery(
            [ [PauliI; PauliI; PauliI; PauliI; PauliI]; 
            [PauliI; PauliX; PauliI; PauliI; PauliI]; 
            [PauliI; PauliI; PauliI; PauliI; PauliZ]; 
            [PauliI; PauliI; PauliX; PauliI; PauliI]; 
            [PauliI; PauliI; PauliZ; PauliI; PauliI]; 
            [PauliZ; PauliI; PauliI; PauliI; PauliI]; 
            [PauliI; PauliI; PauliI; PauliX; PauliI]; 
            [PauliI; PauliI; PauliY; PauliI; PauliI]; 
            [PauliX; PauliI; PauliI; PauliI; PauliI]; 
            [PauliI; PauliI; PauliI; PauliZ; PauliI]; 
            [PauliI; PauliZ; PauliI; PauliI; PauliI]; 
            [PauliI; PauliY; PauliI; PauliI; PauliI]; 
            [PauliI; PauliI; PauliI; PauliI; PauliX]; 
            [PauliY; PauliI; PauliI; PauliI; PauliI]; 
            [PauliI; PauliI; PauliI; PauliI; PauliY]; 
            [PauliI; PauliI; PauliI; PauliY; PauliI] ]
        );
    }

    /// # Summary
    /// Encodes into the ⟦5, 1, 3⟧ quantum code. 
    ///
    /// # Input
    /// ## physRegister
    /// An array of qubits representing an unencoded state.</param>
    /// ## auxQubits
    /// A register of auxillary qubits that will be used to represent the encoded state.</param>
    /// # See Also
    /// - @"microsoft.quantum.canon.FiveQubitCodeDecoder"
    ///
    /// # Remarks
    /// This code was found independently in the following two papers:
    /// C. H. Bennett, D. DiVincenzo, J. A. Smolin and W. K. Wootters, "Mixed state entanglement and quantum error correction,"
    /// Phys. Rev. A, 54 (1996) pp. 3824-3851; https://arxiv.org/abs/quant-ph/9604024 and
    /// R. La?amme, C. Miquel, J. P. Paz and W. H. Zurek, "Perfect quantum error correction code," 
    /// Phys. Rev. Lett. 77 (1996) pp. 198-201; https://arxiv.org/abs/quant-ph/9602019
    operation FiveQubitCodeEncoder(physRegister : Qubit[], auxQubits : Qubit[])  : LogicalRegister
    {
        body {
            FiveQubitCodeEncoderImpl(physRegister, auxQubits);
            
            let logicalRegister = LogicalRegister(physRegister + auxQubits);
            return logicalRegister;
        }
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.FiveQubitCodeEndocer"
    operation FiveQubitCodeDecoder( logicalRegister : LogicalRegister)  : (Qubit[], Qubit[])
    {
        body {
            let physRegister = [logicalRegister[0]];
            let auxQubits = logicalRegister[1..4];

            (Adjoint FiveQubitCodeEncoderImpl)(physRegister, auxQubits);

            return (physRegister, auxQubits);
        }
    }

    /// # Summary
    /// Returns a QECC value representing the ⟦5, 1, 3⟧ code encoder and
    /// decoder with in-place syndrome measurement.
    operation  FiveQubitCode()  : QECC
    {
        body {
            let e = EncodeOp(FiveQubitCodeEncoder);
            let d = DecodeOp(FiveQubitCodeDecoder);
            let s = SyndromeMeasOp(MeasureStabilizerGenerators(
                        [ [ PauliX; PauliZ; PauliZ; PauliX; PauliI ]; 
                        [ PauliI; PauliX; PauliZ; PauliZ; PauliX ];
                        [ PauliX; PauliI; PauliX; PauliZ; PauliZ ];
                        [ PauliZ; PauliX; PauliI; PauliX; PauliZ ] ],
                        _, MeasureWithScratch)
                    );
            let code = QECC(e, d, s);
            return code;
        }
    }

}
