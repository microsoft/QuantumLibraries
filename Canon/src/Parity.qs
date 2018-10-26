// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// This computes the parity of qubits $\ket{q_0} \ket{q_1} \cdots$
    /// in-place, following the pattern
    /// $\ket{q_0} \ket{q_0 \oplus q_1} \ket{q_0 \oplus q_1 \oplus q_2} \cdots$.
    ///
    /// # Input
    /// ## qubits
    /// Array of qubits on which parity is computed 
    /// and stored.
    operation CNOTChain(qubits: Qubit[]) : Unit {
        body (...) {
            for (idxQubit in 0..Length(qubits) - 2){
                CNOT(qubits[idxQubit], qubits[idxQubit + 1]);
            }                    
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// This computes the parity of qubits initially in the state
    /// $\ket{q_0} \ket{q_1} \cdots \ket{q_{\text{target}}}$ in-place,
    /// such that the final state is given by
    /// $\ket{q_0} \ket{q_1 \oplus q_0} \cdots \ket{q_{n - 1} \oplus \cdots \oplus q_0 \oplus q_{\text{target}}}$.
    ///
    /// # Input
    /// ## qubits
    /// Array of qubits on which the parity is computed.
    /// ## targetQubit
    /// Final qubit into which the parity of 'qubits' is XORed.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// CNOTChainTarget(Most(qs), Last(qs));
    /// CNOTChain(qs);
    /// ```
    operation CNOTChainTarget(qubits: Qubit[], targetQubit: Qubit) : Unit {
        body (...) {
            let allQubits = qubits + [targetQubit];
            CNOTChain(allQubits);                 
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// This computes the exclusive-OR of two bits.
    function XOR(bit1: Bool, bit2: Bool) : Bool {
        return (bit1 || bit2) && (not bit1 || not bit2);
    }


}
