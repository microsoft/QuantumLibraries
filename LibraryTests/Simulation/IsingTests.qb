// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Ising {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Canon;
    //open Microsoft.Quantum.Samples.Ising;


    /// Test Ising model anti-ferromagnetic simulation by ZZ correlation function
    operation Ising1DAntiFerromagneticTest() : () {
        body {
            let nSites = 5;
            let adiabaticTime = 100.1;
            let trotterOrder = 1;
            let scheduleSteps = 100;
            let trotterStepSize = adiabaticTime / Float(scheduleSteps);
            let hXamplitude = 1.123;
            let jCamplitude = -0.985;

            // Probabilities obtained from independent simulation
            let probX = [0.498979; 0.49967; 0.499805; 0.49967; 0.498979];
            let probZZ = [0.0000442226; 0.0000213399; 0.0000213399; 0.0000442226];
            using (qubits = Qubit[nSites]) {
                Ising1DStatePrep(qubits);
                (IsingUniformAdiabaticEvolution(nSites, hXamplitude, jCamplitude, adiabaticTime, trotterStepSize, trotterOrder))(qubits);

                for (idxQubit in 0..4) {
                    AssertProb([PauliX], [qubits[idxQubit]], One, probX[idxQubit], "IsingUniformAdiabaticEvolution Qubit X expectation incorrect",  1e-3);
                }
                for (idxQubit in 0..3) {
                    AssertProb([PauliZ; PauliZ], qubits[idxQubit..idxQubit+1], Zero, probZZ[idxQubit], "IsingUniformAdiabaticEvolution Qubit ZZ expectation incorrect",  1e-9);
                }

                ResetAll(qubits);
                Ising1DStatePrep(qubits);
                let hXfinal = Float(0);
                (IsingAdiabaticEvolution_2(nSites, hXamplitude, hXfinal, jCamplitude, adiabaticTime, trotterStepSize, trotterOrder))(qubits);
                for (idxQubit in 0..4) {
                    AssertProb([PauliX], [qubits[idxQubit]], One, probX[idxQubit], "IsingAdiabaticEvolution_2 Qubit X expectation incorrect",  1e-3);
                }
                for (idxQubit in 0..3) {
                    AssertProb([PauliZ; PauliZ], qubits[idxQubit..idxQubit+1], Zero, probZZ[idxQubit], "IsingAdiabaticEvolution_2 Qubit ZZ expectation incorrect",  1e-9);
                }

                ResetAll(qubits);

            }
        }
    }

    // Test Ising model ferromagnetic by ground-state energy estimation at critical point http://dmrg101-tutorial.readthedocs.io/en/latest/tfim.html
    operation Ising1DUniformEstimateTest(): (){
        body {
            let nSites = 6;
            let adiabaticTime = 60.1;
            let trotterOrder = 1;
            let scheduleSteps = 30;
            let qpeStepSize = 0.1;
            let bitsPrecision = 1;

            // FIXME: why are these commented out?
            //let hx = hxLinear(ToDouble(1),_)
            //let jC = GenerateUniform1DJCoupling(nSites,ToDouble(1),_, _)
            //let statePrepUnitary = Ising1DStatePrep
            //let evolutionSchedule = Ising1DSchedule(nSites, hx, jC, trotterOrder)
            ///let qpeUnitary = (evolutionSchedule)(ToDouble(1), qpeStepSize,_)
            ///let qpeUnitary = Ising1DTrotterStepB(nSites, hx,  jC, ToDouble(1), 1, qpeStepSize, _)
            let qpeUnitary = NoOp;
            mutable phaseEst = ToDouble(0);
            using (qubits = Qubit[nSites]) {
                ResetAll(qubits);
                //(AdiabaticStatePrep(adiabaticTime, scheduleSteps, evolutionSchedule, statePrepUnitary))(qubits)
                //set phaseEst = RobustPhaseEstimation(bitsPrecision, OracleToDiscrete(qpeUnitary), qubits) / qpeStepSize      
                ResetAll(qubits);
            }

            //ReadReal(0.5 + phaseEst / ToDouble(100))
            //let (phaseEst, results) = Ising1DUniformEstimate(nSites, hx, jC, adiabaticTime, bitsPrecision, trotterOrder, scheduleSteps, qpeStepSize)

            //ReadReal(0.5 + phaseEst / 100.001)
        }
    }
}
