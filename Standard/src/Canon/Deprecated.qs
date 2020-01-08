// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Logical;

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

    /// # Deprecated
    /// Please use @"microsoft.quantum.logical.xor".
    @Deprecated("Microsoft.Quantum.Logical.Xor")
    function XOR(bit1 : Bool, bit2 : Bool) : Bool {
        return Xor(bit1, bit2);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.applycnotchain".
    @Deprecated("Microsoft.Quantum.Canon.ApplyCNOTChain")
    operation CNOTChain(qubits : Qubit[]) : Unit is Adj + Ctl {
        ApplyCNOTChain(qubits);
    }

}
