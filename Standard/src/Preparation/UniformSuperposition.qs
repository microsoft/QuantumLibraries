// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Creates a uniform superposition over all integers less than a given
    /// index.
    ///
    /// # Description
    /// Given a register of qubits representing a little-endian encoded integer
    /// and initially in the state $\ket{0}$, this operation prepares the
    /// register in the uniform superposition
    /// $$
    /// \begin{align}
    ///    \frac{1}{\sqrt{M}} \sum_{j = 0}^{M - 1} \ket{j}
    /// \end{align}
    /// $$
    /// for a given integer $M$.
    ///
    /// # Input
    /// ## nIndices
    /// The desired number of states $M$ in the uniform superposition.
    /// ## indexRegister
    /// The qubit register that stores the number states in `LittleEndian` format.
    /// This register must be able to store the number $M - 1$, and is assumed to be
    /// initialized in the state $\ket{0}$ (represented by the computational
    /// basis state $\ket{0\cdots 0}).
    ///
    /// # Example
    /// The following example prepares the state $\frac{1}{\sqrt{6}}\sum_{j=0}^{5}\ket{j}$
    /// on $3$ qubits:
    /// ``` Q#
    /// let nIndices = 6;
    /// using (indexRegister = Qubit[3]) {
    ///     PrepareUniformSuperposition(nIndices, LittleEndian(indexRegister));
    ///     // ...
    /// }
    /// ```
    operation PrepareUniformSuperposition(nIndices: Int, indexRegister: LittleEndian) : Unit is Adj + Ctl {
        if(nIndices == 0) {
            fail $"Cannot prepare uniform superposition over {nIndices} state.";
        } elif (nIndices == 1) {
            // Superposition over one state, so do nothing.
        } elif (nIndices == 2) {
            H(indexRegister![0]);
        } else {
            let nQubits = Ceiling(Lg(IntAsDouble(nIndices)));
            if (nQubits > Length(indexRegister!)) {
                fail $"Cannot prepare uniform superposition over {nIndices} states as it is larger than the qubit register.";
            }

            using (flagQubits = Qubit[3]) {
                let targetQubits = indexRegister![0..nQubits - 1];
                let qubits = flagQubits + targetQubits;
                let stateOracle = StateOracle(_PrepareUniformSuperposition(nIndices, nQubits, _, _));

                (StandardAmplitudeAmplification(1, stateOracle, 0))(qubits);

                ApplyToEachCA(X, flagQubits);
            }
        }
    }

    /// # Summary
    /// Implementation step of <xref:microsoft.quantum.canon.prepareuniformsuperposition>
    operation _PrepareUniformSuperposition(nIndices : Int, nQubits : Int, idxFlag : Int, qubits : Qubit[])
    : Unit is Adj + Ctl {
        let targetQubits = qubits[3..3 + nQubits - 1];
        let flagQubit = qubits[0];
        let auxillaryQubits = qubits[1..2];
        let theta = ArcSin(
            Sqrt(IntAsDouble(2 ^ nQubits) / IntAsDouble(nIndices)) *
            Sin(PI() / 6.0)
        );

        ApplyToEachCA(H, targetQubits);
        using (compareQubits = Qubit[nQubits]) {
            within {
                ApplyXorInPlace(nIndices - 1, LittleEndian(compareQubits));
            } apply {
                CompareUsingRippleCarry(LittleEndian(targetQubits), LittleEndian(compareQubits), auxillaryQubits[0]);
                X(auxillaryQubits[0]);
            }
        }
        Exp([PauliY], -theta, [auxillaryQubits[1]]);
        (Controlled X)(auxillaryQubits, flagQubit);
    }

}
