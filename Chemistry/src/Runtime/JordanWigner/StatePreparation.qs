// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;
    

    //newtype JordanWignerInputState = ((Double, Double), Int[]);
    operation PrepareTrialState (stateData : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit {
        let (stateType, terms) = stateData;

        // State type indexing from FermionHamiltonianStatePrep
        // public enum StateType
        //{
        //    Default = 0, Single_Configurational = 1, Sparse_Multi_Configurational = 2, Unitary_Coupled_Cluster = 3
        //}
		if(stateType == 2){
			if (Length(terms) == 0) {
            // Do nothing
			}
			elif (Length(terms) == 1) {
				let (complex, qubitIndices) = terms[0]!;
				PrepareTrialStateSingleSiteOccupation(qubitIndices, qubits);
			}
			else {
				PrepareTrialStateSparseMultiConfigurational(NoOp<Qubit[]>, terms, qubits);
			}
		}
		elif(stateType == 3){
			let nTerms = Length(terms);
			let trotterStepSize = 1.0;

			// The last term is the reference state.
			let referenceState = PrepareTrialState((2, [terms[nTerms-1]]), _);
			
			PrepareTrialStateUnitaryCoupledCluster(referenceState, terms[0..nTerms-2], trotterStepSize, qubits);
		}
    }
    
    
    /// # Summary
    /// Simple state preparation of trial state by occupying
    /// spin-orbitals
    ///
    /// # Input
    /// ## qubitIndices
    /// Indices of qubits to be occupied by electrons.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation PrepareTrialStateSingleSiteOccupation (qubitIndices : Int[], qubits : Qubit[]) : Unit {
        
        body (...) {
            ApplyToEachCA(X, Subarray(qubitIndices, qubits));
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _PrepareTrialStateSingleSiteOccupation_ (qubitIndices : Int[]) : (Qubit[] => Unit : Adjoint, Controlled) {
        
        return PrepareTrialStateSingleSiteOccupation(qubitIndices, _);
    }
    
    
    /// # Summary
    /// Sparse multi-configurational state preparation of trial state by adding excitations
    /// to initial trial state.
    ///
    /// # Input
    /// ## initialStatePreparation
    /// Unitary to prepare initial trial state.
    /// ## excitations
    /// Excitations of initial trial state represented by
    /// the amplitude of the excitation and the qubit indices
    /// the excitation acts on.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation PrepareTrialStateSparseMultiConfigurational (initialStatePreparation : (Qubit[] => Unit), excitations : JordanWignerInputState[], qubits : Qubit[]) : Unit {
        
        let nExcitations = Length(excitations);
        
        //FIXME compile error let coefficientsSqrtAbs = Mapped(Compose(Compose(Sqrt, Fst),Fst), excitations);
        mutable coefficientsSqrtAbs = new Double[nExcitations];
        mutable coefficientsNewComplexPolar = new ComplexPolar[nExcitations];
        mutable applyFlips = new Int[][nExcitations];
        
        for (idx in 0 .. nExcitations - 1) {
            let (x, excitation) = excitations[idx]!;
            set coefficientsSqrtAbs[idx] = Sqrt(AbsComplexPolar(ComplexToComplexPolar(Complex(x))));
            set coefficientsNewComplexPolar[idx] = ComplexPolar(coefficientsSqrtAbs[idx], ArgComplexPolar(ComplexToComplexPolar(Complex(x))));
            set applyFlips[idx] = excitation;
        }
        
        let nBitsIndices = Ceiling(Lg(ToDouble(nExcitations)));
        
        repeat {
            mutable success = false;
            
            using (auxillary = Qubit[nBitsIndices + 1]) {
                using (flag = Qubit[1]) {
                    let multiplexer = MultiplexerBruteForceFromGenerator(nExcitations, LookupFunction(Mapped(_PrepareTrialStateSingleSiteOccupation_, applyFlips)));
                    (StatePreparationComplexCoefficients(coefficientsNewComplexPolar))(BigEndian(auxillary));
                    multiplexer(BigEndian(auxillary), qubits);
                    (Adjoint (StatePreparationPositiveCoefficients(coefficientsSqrtAbs)))(BigEndian(auxillary));
                    (ControlledOnInt(0, X))(auxillary, flag[0]);
                    
                    // if measurement outcome one we prepared required state
                    let outcome = Measure([PauliZ], flag);
                    set success = outcome == One;
                    ResetAll(auxillary);
                    ResetAll(flag);
                }
            }
        }
        until (success)
        fixup {
            ResetAll(qubits);
        }
    }
    
	/// # Summary
    /// Unitary coupled-cluster state preparation of trial state 
    ///
    /// # Input
    /// ## initialStatePreparation
    /// Unitary to prepare initial trial state.
    /// ## excitations
    /// Excitations of initial trial state represented by
    /// the amplitude of the excitation and the qubit indices
    /// the excitation acts on.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation PrepareTrialStateUnitaryCoupledCluster (initialStatePreparation : (Qubit[] => Unit), clusterOperator : JordanWignerInputState[], trotterStepSize : Double, qubits : Qubit[]) : Unit {
        body(...){
			let clusterOperatorGeneratorSystem = JordanWignerClusterOperatorGeneratorSystem(clusterOperator);
			let evolutionGenerator = EvolutionGenerator(JordanWignerClusterOperatorEvolutionSet(), clusterOperatorGeneratorSystem);
			let trotterOrder = 1;
			let simulationAlgorithm = (TrotterSimulationAlgorithm(trotterStepSize, trotterOrder))!;
			let oracle = simulationAlgorithm(1.0, evolutionGenerator, _);

			oracle(qubits);
		}
		adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
	}
}


