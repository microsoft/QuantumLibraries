// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
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
    /// Implements a cascade of CCNOT gates controlled on corresponding bits of two
    /// qubit registers, acting on the next qubit of one of the registers.
    /// Starting from the qubits at position 0 in both registers as controls, CCNOT is
    /// applied to the qubit at position 1 of the target register, then controlled by
    /// the qubits at position 1 acting on the qubit at position 2 in the target register,
    /// etc., ending with an action on the target qubit in position `Length(nQubits)-1`.
    ///
    /// # Input
    /// ## register
    /// Qubit register, only used for controls.
    /// ## targets
    /// Qubit register, used for controls and as target.
    ///
    /// # Remarks
    /// The target qubit register must have one qubit more than the other register.
    operation ApplyCCNOTChain(register : Qubit[], targets : Qubit[])
    : Unit is Adj + Ctl {
        let nQubits = Length(targets);

        EqualityFactI(
            nQubits, Length(register) + 1,
            "Target register must have one more qubit."
        );

        ApplyToEachCA(CCNOT, Zip3(register, Most(targets), Rest(targets)));
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
