// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace SystemTests.Molecules {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Chemistry.JordanWigner;   
    
    // We now use the Q# component of the chemistry library to obtain
    // quantum operations that implement real-time evolution by
    // the chemistry Hamiltonian. Below, we consider two examples.
    // - Trotter simulation algorithm
    // - Qubitization simulation algorithm
    
    // These operations are invoked as orcales in the quantum phase estimation
    // algorithm to extract energy estimates of various eigenstate of the
    // Hamiltonian.
    
    // The returned energy estimate is chosen probabilistically, depending on
    // the overlap of the initial trial state. By default, we greedly
    // fill spin-orbitals to minize the diagonal component of the one-electron
    // energies.
    
    //////////////////////////////////////////////////////////////////////////
    // Using Trotterization //////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    // We can now use Canon's phase estimation algorithms to
    // learn the ground state energy using the above simulation.
    operation GetEnergyByTrotterization (qSharpData : JordanWignerEncodingData, nBitsPrecision : Int, trotterStepSize : Double, trotterOrder : Int) : (Double, Double) {
        
        // The data describing the Hamiltonian for all these steps is contained in
        // `qSharpData`
        let (nSpinOrbitals, fermionTermData, statePrepData, energyOffset) = qSharpData!;
        
        // We use a Product formula, also known as `Trotterization` to
        // simulate the Hamiltonian.
        let (nQubits, (rescaleFactor, oracle)) = TrotterStepOracle(qSharpData, trotterStepSize, trotterOrder);
        
        // The operation that creates the trial state is defined below.
        // By default, greedy filling of spin-orbitals is used.
        let statePrep = PrepareTrialState(statePrepData, _);
        
        // We use the Robust Phase Estimation algorithm
        // of Kimmel, Low and Yoder.
        let phaseEstAlgorithm = RobustPhaseEstimation(nBitsPrecision, _, _);
        
        // This runs the quantum algorithm and returns a phase estimate.
        let estPhase = EstimateEnergy(nQubits, statePrep, oracle, phaseEstAlgorithm);
        
        // We obtain the energy estimate by rescaling the phase estimate
        // with the trotterStepSize. We also add the constant energy offset
        // to the estimated energy.
        let estEnergy = estPhase * rescaleFactor + energyOffset;
        
        // We return both the estimated phase, and the estimated energy.
        return (estPhase, estEnergy);
    }
    
}