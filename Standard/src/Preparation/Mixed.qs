// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;

    /// # Summary
	/// Prepares a qubit in the maximally mixed state.
	///
    /// It prepares the given qubit in the $\boldone / 2$ state by applying the depolarizing channel
    /// $$
    /// \begin{align}
    ///     \Omega(\rho) \mathrel{:=} \frac{1}{4} \sum_{\mu \in \{0, 1, 2, 3\}} \sigma\_{\mu} \rho \sigma\_{\mu}^{\dagger},
    /// \end{align}
    /// $$
    /// where $\sigma\_i$ is the $i$th Pauli operator, and where
    /// $\rho$ is a density operator representing a mixed state.
    ///
    /// # Input
    /// ## qubit
    /// A qubit whose state is to be depolarized in the manner
    /// described above.
    ///
    /// # Remarks
    /// The mixed state $\boldone / 2$ describing the result of
    /// applying this operation to a state implicitly describes
    /// an expectation value over random choices made in this operation.
    /// Thus, for any single application, this operation maps pure states
    /// to pure states, but it acts as described in expectation.
    /// In particular, this operation can be used in process tomography
    /// to measure the *non-unital* components of a channel.
    operation PrepareSingleQubitIdentity(qubit : Qubit) : Unit {
        ApplyPauli([RandomSingleQubitPauli()], [qubit]);
    }

    /// # Summary
    /// Given a register, prepares that register in the maximally mixed state.
	///
    /// The register is prepared in the $\boldone / 2^N$ state by applying the 
	/// complete depolarizing
    /// channel to each qubit, where $N$ is the length of the register.
    ///
    /// # Input
    /// ## register
    /// A register whose state is to be depolarized in the manner
    /// described above.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.preparesinglequbitidentity"
    operation PrepareIdentity(register : Qubit[]) : Unit {
        ApplyToEach(PrepareSingleQubitIdentity, register);
    }

    /// # Summary
	/// Prepares a qubit in the +1 (`Zero`) eigenstate of the given Pauli operator.
	/// If the identity operator is given, then the qubit is prepared in the maximally
	/// mixed state.
	///
    /// If the qubit was initially in the $\ket{0}$ state, this operation prepares the
    /// qubit in the $+1$ eigenstate of a given Pauli operator, or, for `PauliI`,
    /// in the maximally mixed state instead (see <xref:microsoft.quantum.canon.preparesinglequbitidentity>).
    ///
    /// If the qubit was in a state other than $\ket{0}$, this operation applies the following gates:
    /// $H$ for `PauliX`, $HS$ for `PauliY`, $I$ for `PauliZ` and
    /// <xref:microsoft.quantum.canon.preparesinglequbitidentity> for `PauliI`.
    ///
    /// # Input
    /// ## basis
    /// A Pauli operator $P$.
    /// ## qubit
    /// A qubit to be prepared.
    operation PrepareQubit (basis : Pauli, qubit : Qubit) : Unit {
        if (basis == PauliI)
        {
            PrepareSingleQubitIdentity(qubit);
        }
        elif (basis == PauliX)
        {
            H(qubit);
        }
        elif (basis == PauliY)
        {
            H(qubit);
            S(qubit);
        }
    }

}
