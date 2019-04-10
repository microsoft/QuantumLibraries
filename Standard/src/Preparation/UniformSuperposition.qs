// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.AmplitudeAmplification;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Math;
    
    /// # Summary
	/// Creates a uniform superposition over states that encode 0 through `nIndices`.
	/// 
    /// That is, this unitary $U$ creates a uniform superposition over all number states
    /// $0$ to $M-1$, given an input state $\ket{0\cdots 0}$. In other words, 
    /// $$
    /// \begin{align}
    /// U\ket{0}=\frac{1}{\sqrt{M}}\sum_{j=0}^{M-1}\ket{j}.
    /// \end{align}
    /// $$.
    ///
    /// # Input
    /// ## nIndices
    /// The desired number of states $M$ in the uniform superposition.
    /// ## indexRegister
    /// The qubit register that stores the number states in `BigEndian` format.
    /// This register must be able to store the number $M-1$, and is assumed to be
    /// initialized in the state $\ket{0\cdots 0}$.
    ///
    /// # Remarks
    /// ## Example
    /// The following example prepares the state $\frac{1}{\sqrt{6}}\sum_{j=0}^{5}\ket{j}$
    /// on $3$ qubits.
    /// ``` Q#
    /// let nIndices = 6;
    /// using(indexRegister = Qubit[3]) {
    ///     PrepareUniformSuperposition(nIndices, BigEndian(indexRegister));
    /// }
    /// ```
    operation PrepareUniformSuperposition(nIndices: Int, indexRegister: BigEndian) : Unit {
        body (...) {
            if(nIndices == 0) {
                fail $"Cannot prepare uniform superposition over {nIndices} state.";
            }
            elif(nIndices == 1) {
                // Superposition over one state, so do nothing.
            }
            elif(nIndices == 2){
                H(indexRegister![Length(indexRegister!) - 1]);
            }
            else{
                let nQubits = Ceiling(Lg(ToDouble(nIndices)));
                if (nQubits > Length(indexRegister!)){
                    fail $"Cannot prepare uniform superposition over {nIndices} states as it is larger than the qubit register.";
                }
            
                using(flagQubit = Qubit[3]){
                    let targetQubits = indexRegister![Length(indexRegister!) - nQubits..Length(indexRegister!) - 1];
                    
                    let qubits = flagQubit + targetQubits;

                    let stateOracle = StateOracle(PrepareUniformSuperposition_(nIndices, nQubits, _, _));

                    (AmpAmpByOracle(1, stateOracle, 0))(qubits);
                    
                    ApplyToEachCA(X, flagQubit);
                }
            }
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Implementation step of <xref:microsoft.quantum.canon.prepareuniformsuperposition>
    operation PrepareUniformSuperposition_(nIndices: Int, nQubits: Int, idxFlag: Int, qubits: Qubit[]) : Unit {
        body (...) {
            let targetQubits = qubits[3..3 + nQubits-1];
            let flagQubit = qubits[0];
            let auxillaryQubits = qubits[1..2];
            let theta = ArcSin(Sqrt(ToDouble(2^nQubits)/ToDouble(nIndices))*Sin(PI() / 6.0));
            //let theta = PI() * 0.5;

            ApplyToEachCA(H, targetQubits);
            using(compareQubits = Qubit[nQubits]){
                InPlaceXorBE(nIndices - 1, BigEndian(compareQubits));
                ApplyRippleCarryComparatorBE(BigEndian(targetQubits), BigEndian(compareQubits), auxillaryQubits[0]);
                X(auxillaryQubits[0]);
                InPlaceXorBE(nIndices - 1, BigEndian(compareQubits));
            }
            Exp([PauliY], -theta, [auxillaryQubits[1]]);
            (Controlled X)(auxillaryQubits, flagQubit);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }
    
}
