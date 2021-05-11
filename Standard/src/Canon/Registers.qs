// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Uses SWAP gates to Reversed the order of the qubits in
    /// a register.
    ///
    /// # Input
    /// ## register
    /// The qubits order of which should be reversed using SWAP gates
    operation SwapReverseRegister (register : Qubit[]) : Unit {
        body (...) {
            let totalQubits = Length(register);
            let halfTotal = totalQubits / 2;

            for i in 0 .. halfTotal - 1 {
                SWAP(register[i], register[(totalQubits - i) - 1]);
            }
        }

        adjoint self;
        controlled distribute;
        controlled adjoint self;
    }

}
