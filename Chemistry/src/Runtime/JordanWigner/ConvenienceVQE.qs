namespace Microsoft.Quantum.Chemistry.JordanWigner.VQE
{		
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Simulation;


    /// # Summary:
    ///     Private wrapper around PrepareTrialState to make it compatible with EstimateFrequencyA by defining a non-matching adjoint. Note that this private wrapper should not be used outside this context.
    ///     EstimateFrequencyA has built-in emulation on the QuantumSimulator which makes its execution much faster
    ///
    /// # Input
    /// ## inputState
    /// The Jordan-Wigner input required for PrepareTrialState to run
    /// ## qubits
    /// A qubit register
    operation _prepareTrialStateWrapper(inputState : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit is Adj
    {
        body (...)
        {
            Microsoft.Quantum.Diagnostics.AssertAllZero(qubits);
            PrepareTrialState(inputState, qubits);
        }

        adjoint (...)
        {
            ResetAll(qubits);
        }
    }


    /// # Summary:
    ///     Estimate the energy of the molecule by summing the expectation value of individual JW terms.
    ///     Implicitly rely on the spin up-down indexing convention.
    ///
    /// # Input
    /// ## jwHamiltonian
    /// The Jordan-Wigner Hamiltonian object
    /// ## nSamples
    /// The number of samples to use for the estimation of the term expectations
    ///
    /// # Output
    /// The estimated energy of the molecule (does not include constant term)
    operation EstimateEnergy(jwHamiltonian : JordanWignerEncodingData, nSamples : Int) : Double
    {
        // Initialize return values: energy
        mutable energy = 0.;

        // Unpack information and qubit Hamiltonian terms
        let (nQubits, jwTerms, inputState, energyOffset) = jwHamiltonian!;

        // Loop over all qubit Hamiltonian terms
        let generatorsystem = JordanWignerGeneratorSystem(jwTerms);
        let (genInt, genFunc) = generatorsystem!;

        for (idxTerm in 0..genInt-1)
        {
            let term = genFunc(idxTerm);
            let ((idxTermType, coeff), idxFermions) = term!;
            let termType = idxTermType[0];

            let ops = ComputeMeasurementOperators(nQubits, idxFermions, termType);
            let coeffs = ExpandCoefficients(coeff, termType);
            let inputStateUnitary = _prepareTrialStateWrapper(inputState, _);

            let jwTermEnergy = EstimateTermExpectation(inputStateUnitary, ops, coeffs, nQubits, nSamples);
            set energy = energy + jwTermEnergy;
        }

        return energy + energyOffset;
    }
}
