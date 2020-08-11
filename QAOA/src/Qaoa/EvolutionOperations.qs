// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Implements a unitary based on the mixing Hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// ## beta
    /// Vector of coefficents for the unitary based on the mixing Hamiltonian.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.


    operation EvolveWithMixingHamiltonian(qubits: Qubit[], beta: Double) : Unit is Adj + Ctl {
        ApplyToEachCA(R(PauliX, -2.0 * beta, _), qubits);
    }

    /// # Summary
    /// Implements a unitary based on the objective function Hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// ## gamma
    /// Vector of coefficents for the unitary based on the objective function Hamiltonian.
    /// ## oneLocalHamiltonianCoefficients
    /// Array of 1-local coefficents of the objective function Hamiltonian.
    /// ## twoLocalHamiltonianCoefficients
    /// Array of 2-local coefficents of the objective function Hamiltonian.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.

    operation EvolveWithObjectiveHamiltonian(qubits: Qubit[], gamma: Double, oneLocalHamiltonianCoefficients: Double[], twoLocalHamiltonianCoefficients: Double[]) : Unit is Adj + Ctl{
        let numberOfQubits = Length(qubits);
        using (auxiliaryQubit = Qubit()) {
            for(i in 0..numberOfQubits-1) {
                R(PauliZ, 2.0*gamma*oneLocalHamiltonianCoefficients[i],qubits[i]);
            }
            for(i in 0..numberOfQubits-1) {
                for (j in i+1..numberOfQubits-1) {
                    RunPhaseKickback(qubits, auxiliaryQubit, [i,j], 2.0*gamma*twoLocalHamiltonianCoefficients[numberOfQubits*i+j]);
                }
            }
        }
    }
}
