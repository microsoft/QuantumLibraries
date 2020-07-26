namespace Quantum.QAOA {

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
    /// ## h
    /// Array of 1-local coefficents of the objective function Hamiltonian.
    /// ## J
    /// Array of 2-local coefficents of the objective function Hamiltonian.
    /// ## p
    /// Depth of the QAOA circuit.
    ///
    /// # Output
    /// Array of boolean values that represent results of measurements on the QAOA state.
    ///
    /// # References
    /// This implementation in inspired by https://github.com/stephenjordan/qaoa_tsp.

    operation RunQaoa(problemSize: Int, beta: Double[], gamma: Double[], h: Double[], J: Double[], p: Int) : Bool[]
    {
        
        mutable result = new Bool[problemSize];
        using (x = Qubit[problemSize])
        {
            ApplyToEach(H, x);                          // prepare the uniform distribution
            for (i in 0..p-1)
            {
                EvolveWithObjectiveHamiltonian(x, gamma[i], h, J);    // do Exp(-i H_C tz[i])
                EvolveWithMixingHamiltonian(x, beta[i]);            // do Exp(-i H_0 tx[i])
            }
            set result = MeasureAllAndReset(x);                 // measure in the computational basis
        }
        return result;
    }
}
