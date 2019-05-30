// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner.VQE {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Simulation;


    /// # Summary
    /// Wrapper around PrepareTrialState to make it compatible with EstimateFrequencyA by defining an adjoint.
    /// EstimateFrequencyA has built-in emulation feature when targeting the QuantumSimulator which speeds up its execution.
    ///
    /// # Input
    /// ## inputState
    /// The Jordan-Wigner input required for PrepareTrialState to run.
    /// ## qubits
    /// A qubit register.
    operation _prepareTrialStateWrapper(inputState : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit is Adj {

        body (...) {
            Microsoft.Quantum.Diagnostics.AssertAllZero(qubits);
            PrepareTrialState(inputState, qubits);
        }

        // Define a non-matching adjoint body for compliance with EstimateFrequencyA
        adjoint (...) {
            ResetAll(qubits);
        }
    }


    /// # Summary
    /// Estimates the energy of the molecule by summing the energy contributed by the individual Jordan-Wigner terms.
    ///
    /// # Description
    /// This operation implicitly relies on the spin up-down indexing convention.
    ///
    /// # Input
    /// ## jwHamiltonian
    /// The Jordan-Wigner Hamiltonian.
    /// ## nSamples
    /// The number of samples to use for the estimation of the term expectations.
    ///
    /// # Output
    /// The estimated energy of the molecule
    operation EstimateEnergy(jwHamiltonian : JordanWignerEncodingData, nSamples : Int) : Double {

        // Initialize return value
        mutable energy = 0.;

        // Unpack information and qubit Hamiltonian terms
        let (nQubits, jwTerms, inputState, energyOffset) = jwHamiltonian!;

        // Loop over all qubit Hamiltonian terms
	//let dummy = JordanWignerGeneratorSystem(jwTerms);
	let (nTerms, indexFunction) = (JordanWignerGeneratorSystem(jwTerms))!;

        for (idxTerm in 0..nTerms-1) {
        
            let term = indexFunction(idxTerm);
            let ((idxTermType, coeff), idxFermions) = term!;
            let termType = idxTermType[0];

            let ops = MeasurementOperators(nQubits, idxFermions, termType);
            let coeffs = ExpandedCoefficients(coeff, termType);

	    // The private wrapper enables fast emulation during expectation estimation
            let inputStateUnitary = _prepareTrialStateWrapper(inputState, _);

            let jwTermEnergy = EstimateTermExpectation(inputStateUnitary, ops, coeffs, nQubits, nSamples);
            set energy += jwTermEnergy;
        }

        return energy + energyOffset;
    }
}
