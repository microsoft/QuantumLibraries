// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Arithmetic;

    /// # Summary
    /// Represents a particular mixed state that can be prepared on an index
    /// and a garbage register.
    ///
    /// # Input
    /// ## Requirements
    /// Specifies the size of the qubit registers required to prepare the
    /// mixed state represented by this UDT value.
    /// ## Norm
    /// Specifies the 1-norm of the coefficients used to define this mixed
    /// state.
    /// ## Prepare
    /// An operation that, given an index register and a garbage register initially
    /// in the $\ket{0}$ and $\ket{00\cdots 0}$ states (respectively),
    /// prepares the state represented by this UDT value on those registers.
    ///
    /// # See Also
    /// - Microsoft.Quantum.PurifiedMixedState
    newtype MixedStatePreparation = (
        Requirements: MixedStatePreparationRequirements,
        Norm: Double,
        Prepare: ((LittleEndian, Qubit[]) => Unit is Adj + Ctl)
    );

    /// # Summary
    /// Represents the number of qubits required in order to prepare a given
    /// mixed state.
    ///
    /// # Input
    /// ## NTotalQubits
    /// The total number of qubits required by the represented state preparation
    /// operation.
    /// ## NIndexQubits
    /// The number of qubits required for the index register used by the
    /// represented state preparation operation.
    /// ## NGarbageQubits
    /// The number of qubits required for the garbage register used by the
    /// represented state preparation operation.
    ///
    /// # See Also
    /// - Microsoft.Quantum.PurifiedMixedState
    newtype MixedStatePreparationRequirements = (
        NTotalQubits: Int,
        (
            NIndexQubits: Int,
            NGarbageQubits: Int
        )
    );

}
