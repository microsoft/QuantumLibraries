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
            let inputStateUnitary = PrepareTrialStateWrapper(inputState, _);

            let jwTermEnergy = EstimateTermExpectation(inputStateUnitary, ops, coeffs, nQubits, nSamples);
            set energy = energy + jwTermEnergy;
        }

        return energy + energyOffset;
    }
}
