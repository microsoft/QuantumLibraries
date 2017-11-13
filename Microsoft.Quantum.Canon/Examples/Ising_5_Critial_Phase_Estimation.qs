// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Ising {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    /// We may simulate time-evolution to the critical point of this Ising model this way
    /// At the critical point, hXfinal = jZFinal
    function IsingCritialAdiabaticEvolution(nSites: Int, jZFinal: Double, adiabaticTime: Double, trotterStepSize: Double, trotterOrder: Int) : (Qubit[] => () : Adjoint, Controlled) {
        let hXfinal = 2.0 * jZFinal;
        let hXInitial = Float(1);
        return IsingAdiabaticEvolution_2(nSites, hXInitial, hXfinal, jZFinal, adiabaticTime, trotterStepSize, trotterOrder);
    }

    /// At the critical point, we perform phase estimation to exact the ground state energy.
    /// First, we define the unitary on which phase esitmation is performed
    operation IsingCriticalQPEUnitary(nSites: Int, jZFinal: Double, qpeStepSize: Double, qubits: Qubit[]) : () {
        body {
            let hXInitial = Float(1);
            let hXfinal = 2.0 * jZFinal;
            let schedule = Float(1);
            let trotterOrder = 1;
            let simulationAlgorithm = TrotterSimulationAlgorithm(qpeStepSize, trotterOrder);
            let evolutionSet = PauliEvolutionSet();
            let evolutionGenerator = EvolutionGenerator(evolutionSet , IsingEvolutionScheduleImpl(nSites, hXInitial, hXfinal, jZFinal, schedule) );
            simulationAlgorithm(qpeStepSize, evolutionGenerator, qubits);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }


    
    /// Let us define the phase estimation algorithm that we will use
	///     We use the Robust Phase Estimation algorithm 
	///     of Kimmel, Low and Yoder.
    operation IsingCriticalEstimateEnergy(nSites: Int, jZFinal: Double, adiabaticTime: Double, trotterStepSize: Double, trotterOrder: Int, qpeStepSize: Double, nBitsPrecision: Int) : (Double, Result[]) {
        body {
            let qpeOracle = OracleToDiscrete (IsingCriticalQPEUnitary(nSites, jZFinal, qpeStepSize, _) );
            let qpeAlgorithm = RobustPhaseEstimation(nBitsPrecision, _, _);
            let adiabaticEvolution = IsingCritialAdiabaticEvolution(nSites, jZFinal, adiabaticTime, trotterStepSize, trotterOrder);

            mutable phaseEst = Float(0);
            mutable results = new Result[nSites];

            using (qubits = Qubit[nSites]) {
                Ising1DStatePrep(qubits);
                adiabaticEvolution(qubits);
                set phaseEst = qpeAlgorithm(qpeOracle, qubits);
                set results = MultiM(qubits);
                ResetAll(qubits);
            }
            return (phaseEst, results);
        }
    }

    /// Alternatively, we may use the built-in function AdiabaticStateEnergyEstimate
    operation IsingCriticalEstimateEnergy_Builtin(nSites: Int, jZFinal: Double, adiabaticTime: Double, trotterStepSize: Double, trotterOrder: Int, qpeStepSize: Double, nBitsPrecision: Int) : Double {
        body {
            let statePrepUnitary = Ising1DStatePrep;
            let adiabaticUnitary = IsingCritialAdiabaticEvolution(nSites, jZFinal, adiabaticTime, trotterStepSize, trotterOrder);
            let qpeUnitary = IsingCriticalQPEUnitary(nSites, jZFinal, qpeStepSize, _);
            let phaseEstAlgorithm = RobustPhaseEstimation(nBitsPrecision, _, _);

            let phaseEst = AdiabaticStateEnergyEstimate(nSites, statePrepUnitary, adiabaticUnitary, qpeUnitary, phaseEstAlgorithm);
            return phaseEst;
        }
    }

}
