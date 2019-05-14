// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Warnings;
    open Microsoft.Quantum.Arithmetic;

    /// # Deprecated
    /// This operation has been removed.
    operation ApplyQuantumFourierTransformBE(qs : BigEndian) : Unit is Adj + Ctl {
        _Removed(
            "Microsoft.Quantum.Canon.ApplyQuantumFourierTransformBE",
            "ApplyQuantumFourierTransform(BigEndianAsLittleEndian(qs))"
        );
        ApplyQuantumFourierTransform(BigEndianAsLittleEndian(qs));
    }

    operation ApplyQuantumFourierTransformLE(qs : LittleEndian) : Unit is Adj + Ctl {
        _Renamed(
            "Microsoft.Quantum.Canon.ApplyQuantumFourierTransformLE",
            "Microsoft.Quantum.Canon.ApplyQuantumFourierTransform"
        );
        ApplyQuantumFourierTransform(qs);
    }

}
