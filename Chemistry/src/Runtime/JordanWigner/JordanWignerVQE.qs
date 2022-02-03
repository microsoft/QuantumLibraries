// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner.VQE {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Math;


    /// # Summary
    /// Computes the energy associated to a given Jordan-Wigner Hamiltonian term
    ///
    /// # Description
    /// This operation estimates the expectation value associated to each measurement operator and
    /// multiplies it by the corresponding coefficient, using sampling.
    /// The results are aggregated into a variable containing the energy of the Jordan-Wigner term.
    ///
    /// # Input
    /// ## inputStateUnitary
    /// The unitary used for state preparation.
    /// ## ops
    /// The measurement operators of the Jordan-Wigner term.
    /// ## coeffs
    /// The coefficients of the Jordan-Wigner term.
    /// ## nQubits
    /// The number of qubits required to simulate the molecular system.
    /// ## nSamples
    /// The number of samples to use for the estimation of the term expectation.
    ///
    /// # Output
    /// The energy associated to the Jordan-Wigner term.
    operation EstimateTermExpectation(inputStateUnitary : (Qubit[] => Unit is Adj), ops : Pauli[][], coeffs : Double[], nQubits : Int,  nSamples : Int) : Double {
        mutable jwTermEnergy = 0.;

        for (coeff, op) in Zipped(coeffs, ops) {
            // Only perform computation if the coefficient is significant enough
            if (AbsD(coeff) >= 1e-10) {
                // Compute expectation value using the fast frequency estimator, add contribution to Jordan-Wigner term energy
                let termExpectation = EstimateFrequencyA(inputStateUnitary, Measure(op, _), nQubits, nSamples);
                set jwTermEnergy += (2. * termExpectation - 1.) * coeff;
            }
        }

        return jwTermEnergy;
    }


    /// # Summary
    /// Computes all the measurement operators required to compute the expectation of a Jordan-Wigner term.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits required to simulate the molecular system.
    /// ## indices
    /// An array containing the indices of the qubit each Pauli operator is applied to.
    /// ## termType
    /// The type of the Jordan-Wigner term.
    ///
    /// # Output
    /// An array of measurement operators (each being an array of Pauli).
    function MeasurementOperators(nQubits : Int, indices : Int[], termType : Int) : Pauli[][] {
    
        // Compute the size and initialize the array of operators to be returned
        mutable nOps = 0;
        if (termType == 2) {set nOps = 2;}
        elif (termType == 3) {set nOps = 8;}
        else {set nOps = 1;}

        mutable ops = [[], size = nOps];

        // Z and ZZ terms
        if (termType == 0) or (termType == 1) {
            mutable op = ConstantArray(nQubits, PauliI);
            for idx in indices {
                set op w/= idx <- PauliZ;
            }
            set ops w/= 0 <- op;
        }

        // PQRS terms set operators between indices P and Q (resp R and S) to PauliZ
        elif termType == 3 {
            let compactOps = [[PauliX, PauliX, PauliX, PauliX], [PauliY, PauliY, PauliY, PauliY],
                              [PauliX, PauliX, PauliY, PauliY], [PauliY, PauliY, PauliX, PauliX],
                              [PauliX, PauliY, PauliX, PauliY], [PauliY, PauliX, PauliY, PauliX],
                              [PauliY, PauliX, PauliX, PauliY], [PauliX, PauliY, PauliY, PauliX]];

            for iOp in 0..7 {
                mutable compactOp = compactOps[iOp];

                mutable op = ConstantArray(nQubits, PauliI);
                for (idx, pauli) in Zipped(indices, compactOp) {
                    set op w/= idx <- pauli;
                }
                for i in indices[0] + 1..indices[1] - 1 {
                    set op w/= i <- PauliZ;
                }
                for i in indices[2] + 1..indices[3] - 1 {
                    set op w/= i <- PauliZ;
                }
                set ops w/= iOp <- op; 
            }
        }

        // Case of PQ and PQQR terms
        elif termType == 2 {
            let compactOps = [[PauliX, PauliX], [PauliY, PauliY]];

            for iOp in 0..1 {
                mutable compactOp = compactOps[iOp];

                mutable op = ConstantArray(nQubits, PauliI);

                let nIndices = Length(indices);
                set op w/= indices[0] <- compactOp[0];
                set op w/= indices[nIndices-1] <- compactOp[1];
                for i in indices[0] + 1..indices[nIndices - 1] - 1 {
                    set op w/= i <- PauliZ;
                }

                // Case of PQQR term
                if nIndices == 4 {
                     set op w/= indices[1] <- ((indices[0] < indices[1]) and (indices[1] < indices[3])) ? PauliI | PauliZ;
                }
                set ops w/= iOp <- op;
            }
        }

        return ops;
    }


    /// # Summary
    /// Expands the compact representation of the Jordan-Wigner coefficients in order
    /// to obtain a one-to-one mapping between these and Pauli terms.
    ///
    /// # Input
    /// ## coeff
    /// An array of coefficients, as read from the Jordan-Wigner Hamiltonian data structure.
    /// ## termType
    /// The type of the Jordan-Wigner term.
    ///
    /// # Output
    /// Expanded arrays of coefficients, one per Pauli term.
    function ExpandedCoefficients(coeff : Double[], termType : Int) : Double[] {
        // Compute the numbers of coefficients to return
        mutable nCoeffs = 0;
        if (termType == 2) {set nCoeffs = 2;}
        elif (termType == 3) {set nCoeffs = 8;}
        else {set nCoeffs = 1;}

        mutable coeffs = [0.0, size = nCoeffs];

        // Return the expanded array of coefficients
        if (termType == 0) or (termType == 1) {
            set coeffs w/= 0 <- coeff[0];
        }
        elif (termType == 2) or (termType == 3) {
            for i in 0..nCoeffs - 1 {
                set coeffs w/= i <- coeff[i / 2];
            }
        }

        return coeffs;
    }
}
