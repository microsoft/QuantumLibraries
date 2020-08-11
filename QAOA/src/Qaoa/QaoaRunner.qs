// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;


    /// # Summary
    /// Prepares and measures the quantum state in the QAOA.
    ///
    /// # Input
    /// ## problemSize
    /// Number of qubits.
    /// ## beta
    /// Vector of coefficents for the unitary based on the mixing Hamiltonian.
    /// ## gamma
    /// Vector of coefficents for the unitary based on the objective function Hamiltonian.
    /// ## oneLocalHamiltonianCoefficients
    /// Array of 1-local coefficents of the objective function Hamiltonian.
    /// ## twoLocalHamiltonianCoefficients
    /// Array of 2-local coefficents of the objective function Hamiltonian.
    /// ## p
    /// Depth of the QAOA circuit.
    ///
    /// # Output
    /// Array of boolean values that represent results of measurements on the QAOA state.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.
    operation RunQaoa(problemSize: Int, beta: Double[], gamma: Double[], oneLocalHamiltonianCoefficients: Double[], twoLocalHamiltonianCoefficients: Double[], p: Int) : Bool[] {
        
        mutable result = new Bool[problemSize];
        using (qubits = Qubit[problemSize]) {
            ApplyToEach(H, qubits);                         
            for (i in 0..p-1) {
                EvolveWithObjectiveHamiltonian(qubits, gamma[i], oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);
                EvolveWithMixingHamiltonian(qubits, beta[i]);
            }
            set result = MeasureAllAndReset(qubits);                
        }
        return result;
    }
}
