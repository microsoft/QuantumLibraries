namespace Quantum.QAOA {

   open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Implements a unitary based on the mixing hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// ## beta
    /// Vector of coefficents for the unitary based on the mixing Hamiltonian.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.


    operation EvolveWithMixingHamiltonian(qubits: Qubit[], beta: Double) : Unit
    {
        for(i in 0..Length(qubits)-1)
        {
            R(PauliX, -2.0*beta, qubits[i]);
        }
    }

    /// # Summary
    /// Implements a unitary based on the objective function hamiltonian and applies it to qubits.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that will be transformed by a unitary.
    /// ## gamma
    /// Vector of coefficents for the unitary based on the objective function Hamiltonian.
    /// ## h
    /// Array of 1-local coefficents of the objective function Hamiltonian.
    /// ## J
    /// Array of 2-local coefficents of the objective function Hamiltonian.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.

    operation EvolveWithObjectiveHamiltonian(qubits: Qubit[], gamma: Double, h: Double[], J: Double[]) : Unit
    {
	    let numberOfQubits = Length(qubits);
        using (ancillaQubit = Qubit[1])
        {
            for(i in 0..numberOfQubits-1)
            {
                R(PauliZ, 2.0*gamma*h[i],qubits[i]);
            }
            for(i in 0..numberOfQubits-1)
            {
                for (j in i+1..numberOfQubits-1)
                {
                    RunPhaseKickback(qubits, ancillaQubit, [i,j], 2.0*gamma*J[numberOfQubits*i+j]);
                }
            }
        }
    }
}
