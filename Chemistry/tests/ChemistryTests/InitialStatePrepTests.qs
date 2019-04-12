// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Testing;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    
    // Prepare single excitation
    operation PrepareTrialStateCoupledCluster0Test () : Unit {
        
        let nQubits = 6;
        let intTest = [39];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5])];
        
        using (qubits = Qubit[nQubits]) {
            PrepareTrialStateCoupledCluster(NoOp<Qubit[]>, excitations, qubits);
            

            for (idx in IndexRange(excitations)) {
                AssertProbIntBE(intTest[idx], Microsoft.Quantum.Extensions.Math.AbsD(1.0), BigEndian(Reversed(qubits)), 1E-05);
            }
            
            ResetAll(qubits);
        }
    }
    
    
    // Prepare multiple excitations with equal positive weights
    operation PrepareTrialStateCoupledCluster1Test () : Unit {
        
        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5]), JordanWignerInputState((0.1, 0.0), [0, 2, 4]), JordanWignerInputState((0.1, 0.0), [3, 1])];
        
        using (qubits = Qubit[nQubits]) {
            PrepareTrialStateCoupledCluster(NoOp<Qubit[]>, excitations, qubits);
            
            for (idx in IndexRange(excitations)) {
                AssertProbIntBE(intTest[idx], expectedProb[idx], BigEndian(Reversed(qubits)), 1E-05);
            }
            
            ResetAll(qubits);
        }
    }
    
    
    // Prepare multiple excitations with unequal positive weights
    operation PrepareTrialStateCoupledCluster2Test () : Unit {
        
        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [0.047619, 0.190476, 0.761905];
        let excitations = [JordanWignerInputState((0.1, 0.0), [0, 1, 2, 5]), JordanWignerInputState((0.2, 0.0), [0, 2, 4]), JordanWignerInputState((0.4, 0.0), [3, 1])];
        
        using (qubits = Qubit[nQubits]) {
            PrepareTrialStateCoupledCluster(NoOp<Qubit[]>, excitations, qubits);
            
            for (idx in IndexRange(excitations)) {
                AssertProbIntBE(intTest[idx], expectedProb[idx], BigEndian(Reversed(qubits)), 1E-05);
            }
            
            ResetAll(qubits);
        }
    }
    
    
    // Prepare multiple excitations with unequal complex weights
    operation PrepareTrialStateCoupledCluster3Test () : Unit {
        
        let nQubits = 6;
        let intTest = [39, 21, 10];
        let expectedProb = [0.047619, 0.190476, 0.761905];
        let p = [0.3, 0.9, 2.4];
        let excitations = [
            JordanWignerInputState(
                (0.1 * Microsoft.Quantum.Extensions.Math.Cos(p[0]),
                0.1 * Microsoft.Quantum.Extensions.Math.Sin(p[0])),
                [0, 1, 2, 5]
            ),
            
            JordanWignerInputState(
                (0.2 * Microsoft.Quantum.Extensions.Math.Cos(p[1]),
                0.2 * Microsoft.Quantum.Extensions.Math.Sin(p[1])),
                [0, 2, 4]
            ),
            
            JordanWignerInputState(
                (0.4 * Microsoft.Quantum.Extensions.Math.Cos(p[2]), 0.4 * Microsoft.Quantum.Extensions.Math.Sin(p[2])),
                [3, 1]
            )];
        
        using (qubits = Qubit[nQubits]) {
            PrepareTrialStateCoupledCluster(NoOp<Qubit[]>, excitations, qubits);
            
            for (idx in IndexRange(excitations)) {
                AssertProbIntBE(intTest[idx], expectedProb[idx], BigEndian(Reversed(qubits)), 1E-05);
            }
            
            ResetAll(qubits);
        }
    }
    
    
    // Prepare multiple excitations with complex weights
    operation PrepareTrialStateCoupledCluster4Test () : Unit {
        
        let nQubits = 1;
        let intTest = [39, 21, 10];
        let phase = 2.453;
        let excitations = [JordanWignerInputState((0.1, 0.0), new Int[0]), JordanWignerInputState((0.1 * Microsoft.Quantum.Extensions.Math.Cos(phase), 0.1 * Microsoft.Quantum.Extensions.Math.Sin(phase)), [0])];
        
        using (qubits = Qubit[nQubits]) {
            PrepareTrialStateCoupledCluster(NoOp<Qubit[]>, excitations, qubits);
            AssertPhase(-phase / 2.0, qubits[0], 1E-09);
            ResetAll(qubits);
        }
    }
    
}


