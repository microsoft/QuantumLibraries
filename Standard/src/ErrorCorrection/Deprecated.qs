// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Canon;

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintobitflipcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.EncodeIntoBitFlipCode")
    operation BitFlipEncoder(physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        return Microsoft.Quantum.ErrorCorrection.EncodeIntoBitFlipCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintobitflipcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.DecodeFromBitFlipCode")
    operation BitFlipDecoder(logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        return Microsoft.Quantum.ErrorCorrection.DecodeFromBitFlipCode(logicalRegister);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintofivequbitcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.EncodeIntoFiveQubitCode")
    operation FiveQubitCodeEncoder(physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        return EncodeIntoFiveQubitCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.DecodeFromFiveQubitCode")
    operation FiveQubitCodeDecoder(logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        return DecodeFromFiveQubitCode(logicalRegister);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.EncodeIntoSteaneCode")
    operation SteaneCodeEncoder (physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        return EncodeIntoSteaneCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    @Deprecated("Microsoft.Quantum.ErrorCorrection.DecodeFromSteaneCode")
    operation SteaneCodeDecoder (logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        return DecodeFromSteaneCode(logicalRegister);
    }

}
