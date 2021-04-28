// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Measures the identity operator on a register
    /// of qubits.
    ///
    /// # Input
    /// ## register
    /// The register to be measured.
    ///
    /// # Output
    /// The result value `Zero`.
    ///
    /// # Remarks
    /// Since $\boldone$ has only the eigenvalue $1$, and does not
    /// have a negative eigenvalue, this operation always returns
    /// `Zero`, corresponding to the eigenvalue $+1 = (-1)^0$,
    /// and does not cause a collapse of the state of `register`.
    ///
    /// On its own, this operation is not useful, but is helpful
    /// in the context of process tomography, as it provides
    /// information about the trace preservation of a quantum process.
    /// In particular, a target machine with lossy measurement should
    /// replace this operation by an actual measurement of $\boldone$.
    operation MeasureIdentity (register : Qubit[]) : Result {
        return Zero;
    }

    /// # Summary
    /// Performs a single-qubit process tomography measurement in the Pauli
    /// basis, given a particular channel of interest.
    ///
    /// # Input
    /// ## preparation
    /// The Pauli basis element $P$ in which a qubit is to be prepared.
    /// ## measurement
    /// The Pauli basis element $Q$ in which a qubit is to be measured.
    /// ## channel
    /// A single qubit channel $\Lambda$ whose behavior is being estimated
    /// with process tomography.
    ///
    /// # Output
    /// The Result `Zero` with probability
    /// $$
    /// \begin{align}
    ///     \Pr(\texttt{Zero} | \Lambda; P, Q) = \operatorname{Tr}\left(
    ///         \frac{\boldone + Q}{2} \Lambda\left[
    ///             \frac{\boldone + P}{2}
    ///         \right]
    ///     \right).
    /// \end{align}
    /// $$
    ///
    /// # Remarks
    /// The distribution over results returned by this operation is a special
    /// case of two-qubit state tomography. Let $\rho = J(\Lambda) / 2$ be
    /// the Choi–Jamiłkowski state for $\Lambda$. Then, the distribution above
    /// is identical to
    /// $$
    /// \begin{align}
    ///     \Pr(\texttt{Zero} | \rho; M) = \operatorname{Tr}(M \rho),
    /// \end{align}
    /// $$
    /// where $M = 2 (\boldone + P)^\mathrm{T} / 2 \cdot (\boldone + Q) / 2$
    /// is the effective measurement corresponding to $P$ and $Q$.
    operation SingleQubitProcessTomographyMeasurement (preparation : Pauli, measurement : Pauli, channel : (Qubit => Unit)) : Result     {
        mutable result = Zero;

        use qubit = Qubit();
        PreparePauliEigenstate(preparation, qubit);
        channel(qubit);
        set result = Measure([measurement], [qubit]);
        Reset(qubit);

        return result;
    }

}
