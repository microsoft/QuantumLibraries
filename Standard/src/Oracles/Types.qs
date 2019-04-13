// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Oracles {

    /// # Summary
    /// Represents a reflection oracle.
    ///
    /// A reflection oracle, $O$, has inputs:
    /// - The phase $\phi$ by which to rotate the reflected subspace.
    /// - The qubit register on which to perform the given reflection.
    ///
    /// # Remarks
    /// This oracle $O = \boldone - (1 - e^{i \phi}) \ket{\psi}\bra{\psi}$
    /// performs a partial reflection by a phase $\phi$ about a single pure state
    /// $\ket{\psi}$.
    newtype ReflectionOracle = ((Double, Qubit[]) => Unit : Adjoint, Controlled);

    // This oracle O|s>_a|ψ>_s = λ |t>_a U |ψ>_s + ... acts on the ancilla state |s>_a to implement the unitary U on any system state |ψ>_s with amplitude λ in the |t>_a basis.

    /// # Summary
    /// Represents an oracle for oblivious amplitude amplification.
    ///
    /// The inputs to the oracle $O$ are:
    /// - The ancilla register $a$ that $O$ acts on.
    /// - The system register $s$ on which the desired unitary $U$ is applied, post-selected on register $a$ being in state $\ket{t}\_a$.
    ///
    /// # Remarks
    /// This oracle defined by
    /// $$
    ///O\ket{s}\_a\ket{\psi}\_s= \lambda\ket{t}\_a U \ket{\psi}\_s + \sqrt{1-|\lambda|^2}\ket{t^\perp}\_a\cdots
    /// $$
    /// acts on the ancilla state $\ket{s}\_a$ to implement the unitary $U$ on any system state $\ket{\psi}\_s$ with amplitude $\lambda$ in the basis flagged by $\ket{t}\_a$.
    /// The first parameter is the qubit register of $\ket{s}\_a$. The second parameter is the qubit register of $\ket{\psi}\_s$.
    newtype ObliviousOracle = ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled);

    /// # Summary
    /// Represents an oracle for state preparation.
    ///
    /// The inputs to the oracle $O$ are:
    /// - An integer indexing the flag qubit $f$.
    /// - The system register $s$ that will store the desired quantum state $\ket{\psi}\_s$.
    ///
    /// # Remarks
    /// This oracle defined by
    /// $$
    /// O\ket{0}\_{f}\ket{0}\_s= \lambda\ket{1}\_f\ket{\psi}\_s + \sqrt{1-|\lambda|^2}\ket{0}\_f\cdots,
    /// $$
    /// acts on the on computational basis state $\ket{0}\_{f}\ket{0}\_s$ to create the target state $\ket{\psi}\_s$ with amplitude $\lambda$ in the basis flagged by $\ket{1}\_f$.
    /// The first parameter is an index to the qubit register of $\ket{0}\_f$. The second parameter encompassed both registers.
    newtype StateOracle = ((Int, Qubit[]) => Unit : Adjoint, Controlled);

    /// # Summary
    /// Represents an oracle for deterministic state preparation.
    ///
    /// The input to the oracle $O$ is:
    /// - The register that will store the desired quantum state $\ket{\psi}\_s$.
    ///
    /// # Remarks
    /// This oracle defined by $O\ket{0}=\ket{\psi}$ acts on the on computational basis state $\ket{0}$ to create the state $\ket{\psi}$.
    /// The first parameter is the qubit register of $\ket{\psi}$.
    newtype DeterministicStateOracle = (Qubit[] => Unit : Adjoint, Controlled);


    /// # Summary
    /// Represents a discrete-time oracle.
	///
	/// This is an oracle that implements $U^m$ for a fixed operation $U$
    /// and a non-negative integer $m$.
    newtype DiscreteOracle = ((Int, Qubit[]) => Unit : Adjoint, Controlled);

    /// # Summary
    /// Represents a continuous-time oracle.
	///
	/// This is an oracle that implements
    /// $U(\delta t) : \ket{\psi(t)} \mapsto \ket{\psi(t + \delta t)}$
    /// for all times $t$, where $U$ is a fixed operation, and where
    /// $\delta t$ is a non-negative real number.
    newtype ContinuousOracle = ((Double, Qubit[]) => Unit : Adjoint, Controlled);

}


