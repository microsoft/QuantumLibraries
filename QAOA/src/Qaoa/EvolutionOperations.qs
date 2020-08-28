// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Implements a unitary based on the mixing Hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// ## beta
    /// Coefficent for the unitary based on the mixing Hamiltonian.
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.
    operation EvolveWithMixingHamiltonian(beta: Double, qubits: Qubit[]) : Unit is Adj + Ctl {
        ApplyToEachCA(R(PauliX, -2.0 * beta, _), qubits);
    }

    /// # Summary
    /// Implements a unitary based on the objective function Hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## gamma
    /// Coefficent for the unitary based on the objective function Hamiltonian.
    /// ## oneLocalHamiltonianCoefficients
    /// Array of 1-local coefficents of the objective function Hamiltonian. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n. The i-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z.
    /// ## twoLocalHamiltonianCoefficients
    /// Array of 2-local coefficents of the objective function Hamiltonian. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n^2. The (i*n+j)-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z\sigma_j^z (other coefficients in an array can take any value of type double).
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.
    operation EvolveWithObjectiveHamiltonian(gamma: Double, oneLocalHamiltonianCoefficients: Double[], twoLocalHamiltonianCoefficients: Double[], qubits: Qubit[]) : Unit is Adj + Ctl{
        let numberOfQubits = Length(qubits);
        for(i in 0..numberOfQubits-1) {
            R(PauliZ, 2.0*gamma*oneLocalHamiltonianCoefficients[i],qubits[i]);
        }
        for(i in 0..numberOfQubits-1) {
            for (j in i+1..numberOfQubits-1) {
                RunPhaseKickback(2.0*gamma*twoLocalHamiltonianCoefficients[i*numberOfQubits+j], Subarray([i,j], qubits));
            }
        }
    }
}
