// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Simulation;

    //newtype JordanWignerInputState = ((Double, Double), Int[]);
    operation PrepareTrialState (stateData : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit {
        let (stateType, terms) = stateData;

        // State type indexing from FermionHamiltonianStatePrep
        // public enum StateType
        //{
        //    Default = 0, Single_Configurational = 1, Sparse_Multi_Configurational = 2, Unitary_Coupled_Cluster = 3
        //}
        if stateType == 2 {
            if IsEmpty(terms) {
                // Do nothing, as there are no terms to prepare.
            } elif Length(terms) == 1 {
                let (coefficient, qubitIndices) = terms[0]!;
                PrepareSingleConfigurationalStateSingleSiteOccupation(qubitIndices, qubits);
            } else {
                PrepareSparseMultiConfigurationalState(NoOp, terms, qubits);
            }
        } elif stateType == 3 {
            let nTerms = Length(terms);
            let trotterStepSize = 1.0;

            // The last term is the reference state.
            let referenceState = PrepareTrialState((2, [terms[nTerms - 1]]), _);

            PrepareUnitaryCoupledClusterState(referenceState, terms[...nTerms - 2], trotterStepSize, qubits);
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
    operation PrepareSingleConfigurationalStateSingleSiteOccupation (qubitIndices : Int[], qubits : Qubit[])
    : Unit is Adj + Ctl {
        ApplyToEachCA(X, Subarray(qubitIndices, qubits));
    }


    internal function _PrepareSingleConfigurationalStateSingleSiteOccupation(qubitIndices : Int[]) : (Qubit[] => Unit is Adj + Ctl) {
        return PrepareSingleConfigurationalStateSingleSiteOccupation(qubitIndices, _);
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
    operation PrepareSparseMultiConfigurationalState(
        initialStatePreparation : (Qubit[] => Unit),
        excitations : JordanWignerInputState[],
        qubits : Qubit[]
    )
    : Unit {
        let nExcitations = Length(excitations);

        //FIXME compile error let coefficientsSqrtAbs = Mapped(Compose(Compose(Sqrt, Fst),Fst), excitations);
        mutable coefficientsSqrtAbs = new Double[nExcitations];
        mutable coefficientsNewComplexPolar = new ComplexPolar[nExcitations];
        mutable applyFlips = new Int[][nExcitations];

        for idx in 0 .. nExcitations - 1 {
            let (x, excitation) = excitations[idx]!;
            set coefficientsSqrtAbs w/= idx <- Sqrt(AbsComplexPolar(ComplexAsComplexPolar(Complex(x))));
            set coefficientsNewComplexPolar w/= idx <- ComplexPolar(coefficientsSqrtAbs[idx], ArgComplexPolar(ComplexAsComplexPolar(Complex(x))));
            set applyFlips w/= idx <- excitation;
        }

        let nBitsIndices = Ceiling(Lg(IntAsDouble(nExcitations)));

        repeat {
            mutable success = false;
            use auxillary = Qubit[nBitsIndices + 1];
            use flag = Qubit();

            let multiplexer = MultiplexerBruteForceFromGenerator(nExcitations, LookupFunction(Mapped(_PrepareSingleConfigurationalStateSingleSiteOccupation, applyFlips)));
            PrepareArbitraryStateCP(coefficientsNewComplexPolar, LittleEndian(auxillary));
            multiplexer(LittleEndian(auxillary), qubits);
            Adjoint PrepareArbitraryStateD(coefficientsSqrtAbs, LittleEndian(auxillary));
            ControlledOnInt(0, X)(auxillary, flag);

            // if measurement outcome one we prepared required state
            let outcome = M(flag);
            set success = outcome == One;
            ResetAll(auxillary);
            Reset(flag);
        }
        until success
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
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation PrepareUnitaryCoupledClusterState (initialStatePreparation : (Qubit[] => Unit), clusterOperator : JordanWignerInputState[], trotterStepSize : Double, qubits : Qubit[]) : Unit {
        let clusterOperatorGeneratorSystem = JordanWignerClusterOperatorGeneratorSystem(clusterOperator);
        let evolutionGenerator = EvolutionGenerator(JordanWignerClusterOperatorEvolutionSet(), clusterOperatorGeneratorSystem);
        let trotterOrder = 1;
        let simulationAlgorithm = (TrotterSimulationAlgorithm(trotterStepSize, trotterOrder))!;
        let oracle = simulationAlgorithm(1.0, evolutionGenerator, _);
        initialStatePreparation(qubits);
        oracle(qubits);
    }
}
