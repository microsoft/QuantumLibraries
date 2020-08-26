// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;


    /// # Summary
    /// Prepares and measures the quantum state in the QAOA.
    ///
    /// # Input
    /// ## problemSize
    /// Number of qubits.
    /// ## betas
    /// Vector of coefficents for the unitary based on the mixing Hamiltonian.
    /// ## gammas
    /// Vector of coefficents for the unitary based on the objective function Hamiltonian.
    /// ## oneLocalHamiltonianCoefficients
    /// Array of 1-local coefficents of the objective function Hamiltonian.
    /// ## twoLocalHamiltonianCoefficients
    /// Array of 2-local coefficents of the objective function Hamiltonian.
    /// ## p
    /// A parameter related to the depth of a QAOA circuit. It corresponds to the number of times evolution operators are applied.
    ///
    /// # Output
    /// Array of boolean values that represent results of measurements on the QAOA state.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.
    operation RunQaoa(problemSize: Int, betas: Double[], gammas: Double[], oneLocalHamiltonianCoefficients: Double[], twoLocalHamiltonianCoefficients: Double[], p: Int) : Bool[] {
        
        mutable result = new Bool[problemSize];
        using (qubits = Qubit[problemSize]) {
            ApplyToEach(H, qubits);                         
            for ((beta, gamma) in Zip(betas, gammas)) {
                EvolveWithObjectiveHamiltonian(qubits, gamma, oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);
                EvolveWithMixingHamiltonian(qubits, beta);
            }
            return ResultArrayAsBoolArray(ForEach(MResetZ, qubits));               
        }
    }
}
