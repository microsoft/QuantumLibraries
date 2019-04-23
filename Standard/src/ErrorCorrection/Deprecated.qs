// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Canon;

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintobitflipcode".
    operation BitFlipEncoder(physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        _Renamed("Microsoft.Quantum.Canon.BitFlipEncoder", "Microsoft.Quantum.ErrorCorrection.EncodeIntoBitFlipCode");
        return Microsoft.Quantum.ErrorCorrection.EncodeIntoBitFlipCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintobitflipcode".
    operation BitFlipDecoder(logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        _Renamed("Microsoft.Quantum.Canon.BitFlipDecoder", "Microsoft.Quantum.ErrorCorrection.DecodeFromBitFlipCode");
        return Microsoft.Quantum.ErrorCorrection.DecodeFromBitFlipCode(logicalRegister);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.encodeintofivequbitcode".
    operation FiveQubitCodeEncoder(physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        _Renamed("Microsoft.Quantum.Canon.FiveQubitCodeEncoder", "Microsoft.Quantum.ErrorCorrection.EncodeIntoFiveQubitCode");
        return EncodeIntoFiveQubitCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    operation FiveQubitCodeDecoder(logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        _Renamed("Microsoft.Quantum.Canon.FiveQubitCodeDecoder", "Microsoft.Quantum.ErrorCorrection.DecodeFromFiveQubitCode");
        return DecodeFromFiveQubitCode(logicalRegister);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    operation SteaneCodeEncoder (physRegister : Qubit[], auxQubits : Qubit[]) : LogicalRegister {
        _Renamed("Microsoft.Quantum.Canon.SteaneCodeEncoder", "Microsoft.Quantum.ErrorCorrection.EncodeIntoSteaneCode");
        return EncodeIntoSteaneCode(physRegister, auxQubits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.errorcorrection.decodefromfivequbitcode".
    operation SteaneCodeDecoder (logicalRegister : LogicalRegister) : (Qubit[], Qubit[]) {
        _Renamed("Microsoft.Quantum.Canon.SteaneCodeDecoder", "Microsoft.Quantum.ErrorCorrection.DecodeFromSteaneCode");
        return DecodeFromSteaneCode(logicalRegister);
    }

}
