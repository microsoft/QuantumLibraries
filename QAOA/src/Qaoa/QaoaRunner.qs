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
    /// Array of 1-local coefficents of the objective function Hamiltonian. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n. The i-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z.
    /// ## twoLocalHamiltonianCoefficients
    /// Array of 2-local coefficents of the objective function Hamiltonian. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n^2. The (i*n+j)-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z\sigma_j^z (other coefficients in an array can take any value of type double).
    ///
    /// # Output
    /// Array of boolean values that represent results of measurements on the QAOA state.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.
    ///
    /// # Examples
    /// Suppose we want to solve a MAXCUT problem on a simple weighted graph with two connected vertices. Then, problemSize = 2, the corresponding quantum Hamiltonian is \sigma_0^z\sigma_1^z and is expressed as oneLocalHamiltonianCoefficients = new Double[] {0, 0} and twoLocalHamiltonianCoefficients = new Double[]{0, 1, 0, 0}. betas and gammas are arrays of p elements each (their actual values are difficult to choose upfront and p is related to the depth of the QAOA circuit).
    operation RunQaoa(problemSize: Int, betas: Double[], gammas: Double[], oneLocalHamiltonianCoefficients: Double[], twoLocalHamiltonianCoefficients: Double[]) : Bool[] {
        
        mutable result = new Bool[problemSize];
        using (qubits = Qubit[problemSize]) {
            ApplyToEach(H, qubits);                         
            for ((beta, gamma) in Zip(betas, gammas)) {
                EvolveWithObjectiveHamiltonian(gamma, oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients, qubits);
                EvolveWithMixingHamiltonian(beta, qubits);
            }
            return ResultArrayAsBoolArray(ForEach(MResetZ, qubits));               
        }
    }
}
