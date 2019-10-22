// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;

    /// # Deprecated
    /// This operation has been removed.
    @Deprecated("ApplyQuantumFourierTransform(BigEndianAsLittleEndian(qs))")
    operation ApplyQuantumFourierTransformBE(qs : BigEndian) : Unit is Adj + Ctl {
        ApplyQuantumFourierTransform(BigEndianAsLittleEndian(qs));
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.applyquantumfouriertransform".
    @Deprecated("Microsoft.Quantum.Canon.ApplyQuantumFourierTransform")
    operation ApplyQuantumFourierTransformLE(qs : LittleEndian) : Unit is Adj + Ctl {
        ApplyQuantumFourierTransform(qs);
    }

}
