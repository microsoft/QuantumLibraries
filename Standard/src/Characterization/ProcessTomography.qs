// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Preparation;

    /// # Summary
	/// Jointly measures a register of qubits in the Pauli Z basis.
	///
    /// In other words, measures the operation $Z \otimes Z \otimes \cdots \otimes Z$ on
    /// a given register.
    ///
    /// # Input
    /// ## register
    /// The register to be measured.
    ///
    /// # Output
    /// The result of measuring $Z \otimes Z \otimes \cdots \otimes Z$.
    operation MeasureAllZ (register : Qubit[]) : Result {
        return Measure(ConstantArray(Length(register), PauliZ), register);
    }

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
    /// Given a preparation and measurement, estimates the frequency
    /// with which that measurement succeeds (returns `Zero`) by
    /// performing a given number of trials.
    ///
    /// # Input
    /// ## preparation
    /// An operation $P$ that prepares a given state $\rho$ on
    /// its input register.
    /// ## measurement
    /// An operation $M$ representing the measurement of interest.
    /// ## nQubits
    /// The number of qubits on which the preparation and measurement
    /// each act.
    /// ## nMeasurements
    /// The number of times that the measurement should be performed
    /// in order to estimate the frequency of interest.
    ///
    /// # Output
    /// An estimate $\hat{p}$ of the frequency with which
    /// $M(P(\ket{00 \cdots 0}\bra{00 \cdots 0}))$ returns `Zero`,
    /// obtained using the unbiased binomial estimator $\hat{p} =
    /// n\_{\uparrow} / n\_{\text{measurements}}$, where $n\_{\uparrow}$ is
    /// the number of `Zero` results observed.
    ///
    /// This is particularly important on target machines which respect
    /// physical limitations, such that probabilities cannot be measured.
    operation EstimateFrequency (preparation : (Qubit[] => Unit), measurement : (Qubit[] => Result), nQubits : Int, nMeasurements : Int) : Double
    {
        mutable nUp = 0;

        for (idxMeasurement in 0 .. nMeasurements - 1) {
            using (register = Qubit[nQubits]) {
                preparation(register);
                let result = measurement(register);

                if (result == Zero) {
                    // NB!!!!! This reverses Zero and One to use conventions
                    //         common in the QCVV community. That is confusing
                    //         but is confusing with an actual purpose.
                    set nUp = nUp + 1;
                }

                // NB: We absolutely must reset here, since preparation()
                //     and measurement() can each use randomness internally.
                ApplyToEach(Reset, register);
            }
        }

        return ToDouble(nUp) / ToDouble(nMeasurements);
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
    operation SingleQubitProcessTomographyMeasurement (preparation : Pauli, measurement : Pauli, channel : (Qubit => Unit)) : Result
    {
        mutable result = Zero;
        
        using (register = Qubit[1])
        {
            let qubit = register[0];
            PrepareQubit(preparation, qubit);
            channel(qubit);
            set result = Measure([measurement], [qubit]);
            Reset(qubit);
        }
        
        return result;
    }
    
}


