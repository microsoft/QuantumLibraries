// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Computes the parity of a register of qubits in-place.
    ///
    /// # Description
    /// This operation transforms the state of its input as
    /// $$
    /// \begin{align}
    ///     \ket{q_0} \ket{q_1} \cdots \ket{q_{n - 1}} & \mapsto
    ///     \ket{q_0} \ket{q_0 \oplus q_1} \ket{q_0 \oplus q_1 \oplus q_2} \cdots
    ///         \ket{q_0 \oplus \cdots \oplus q_{n - 1}}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## qubits
    /// Array of qubits whose parity is to be computed and stored.
    operation ApplyCNOTChain(qubits : Qubit[]) : Unit is Adj + Ctl {
        ApplyToEachCA(CNOT, Zip(Most(qubits), Rest(qubits)));
    }

    /// # Summary
    /// Computes the parity of an array of qubits into a target qubit.
    ///
    /// # Description
    /// If the array is initially in the state
    /// $\ket{q_0} \ket{q_1} \cdots \ket{q_{\text{target}}}$,
    /// the final state is given by
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
    /// ```qsharp
    /// ApplyCNOTChainWithTarget(Most(qs), Last(qs));
    /// ```
    /// and
    /// ```qsharp
    /// ApplyCNOTChain(qs);
    /// ```
    operation ApplyCNOTChainWithTarget(qubits : Qubit[], targetQubit : Qubit)
    : Unit is Adj + Ctl {
        ApplyCNOTChain(qubits + [targetQubit]);
    }

}
