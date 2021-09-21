// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Logical;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Logical.Xor".
    @Deprecated("Microsoft.Quantum.Logical.Xor")
    function XOR(bit1 : Bool, bit2 : Bool) : Bool {
        return Xor(bit1, bit2);
    }

    @Deprecated("Microsoft.Quantum.Canon.ApplyCNOTChain")
    operation CNOTChain(qubits : Qubit[]) : Unit is Adj + Ctl {
        ApplyCNOTChain(qubits);
    }

    @Deprecated("Microsoft.Quantum.Canon.ApplyCNOTChainWithTarget")
    operation CNOTChainTarget(qubits : Qubit[], targetQubit : Qubit)
    : Unit is Adj + Ctl {
        ApplyCNOTChainWithTarget(qubits, targetQubit);
    }

}
