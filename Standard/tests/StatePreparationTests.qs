// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics as Diag;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Preparation;

    // number of qubits, abs(amplitude), phase
    newtype StatePreparationTestCase = (Int, Double[], Double[]);


    @Diag.Test("QuantumSimulator")
    operation TestStatePreparationPositiveCoefficients() : Unit {

        let tolerance = 1E-09;

        let testCases = [
            // Test positive coefficients
            StatePreparationTestCase(1, [0.773761, 0.633478], [0.0, 0.0]),
            StatePreparationTestCase(2, [0.183017, 0.406973, 0.604925, 0.659502], [0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(3, [0.0986553, 0.359005, 0.465689, 0.467395, 0.419893, 0.118445, 0.461883, 0.149609], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(4, [0.271471, 0.0583654, 0.11639, 0.36112, 0.307383, 0.193371, 0.274151, 0.332542, 0.130172, 0.222546, 0.314879, 0.210704, 0.212429, 0.245518, 0.30666, 0.22773], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),

            // Test negative coefficients; should give same probabilities as positive coefficients
            StatePreparationTestCase(1, [-0.773761, 0.633478], [0.0, 0.0]),
            StatePreparationTestCase(2, [0.183017, -0.406973, 0.604925, 0.659502], [0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(3, [0.0986553, -0.359005, 0.465689, -0.467395, 0.419893, 0.118445, -0.461883, 0.149609], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(4, [-0.271471, 0.0583654, 0.11639, 0.36112, -0.307383, 0.193371, -0.274151, 0.332542, 0.130172, 0.222546, 0.314879, -0.210704, 0.212429, 0.245518, -0.30666, -0.22773], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),

            // Test unnormalized coefficients
            StatePreparationTestCase(3, [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445, 0.461883, 0.149609], new Double[0]),

            // Test missing coefficients
            StatePreparationTestCase(3, [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445], new Double[0])
        ];
        
        // Loop over multiple qubit tests
        for (testCase in testCases) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCase!;
            let nCoefficients = Length(coefficientsAmplitude);
            
            // Test negative coefficients. Should give same results as positive coefficients
            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                let op = StatePreparationPositiveCoefficients(coefficientsAmplitude);
                op(qubitsLE);
                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);
                
                for (idxCoeff in IndexRange(coefficientsAmplitude)) {
                    let amp = normalizedCoefficients[idxCoeff];
                    let prob = amp * amp;
                    AssertProbInt(idxCoeff, prob, qubitsLE, tolerance);
                }
                
                ResetAll(qubits);
            }
        }
    }
    
    
    // Test phase factor on 1-qubit uniform superposition.
    @Diag.Test("QuantumSimulator")
    operation TestStatePreparationComplexCoefficientsQubitPhase() : Unit {
        
        let tolerance = 1E-09;

        let testCases = [
            // Test phase factor on uniform superposition
            StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.01]),
            StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.05])
        ];
        
        // Loop over tests
        for (testCase in testCases) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCase!;
            let nCoefficients = Length(coefficientsAmplitude);
            
            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                let coefficients = Mapped(ComplexPolar, Zip(coefficientsAmplitude, coefficientsPhase));
                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);
                
                // Test phase factor on uniform superposition
                let phase = 0.5 * (coefficientsPhase[0] - coefficientsPhase[1]);
                let amp = normalizedCoefficients[0];
                let prob = amp * amp;
                let op = StatePreparationComplexCoefficients(coefficients);
                op(qubitsLE);
                AssertProbInt(0, prob, qubitsLE, tolerance);
                AssertProbInt(1, prob, qubitsLE, tolerance);
                Diag.AssertPhase(phase, (qubitsLE!)[0], tolerance);
                ResetAll(qubits);
            }
        }
    }
    
    
    // Test probabilities and phases factor of multi-qubit uniform superposition.
    @Diag.Test("QuantumSimulator")
    operation TestStatePreparationComplexCoefficientsMultiQubitPhase() : Unit {
        
        let tolerance = 1E-09;

        let testCases = [
            // Test probability and phases of uniform superposition
            StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.01]),
            StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], ConstantArray(8, PI())),
            StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01]),
            StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445, 0.461883, 0.149609])
        ];
        
        // Loop over tests
        for (testCase in testCases) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCase!;
            let nCoefficients = Length(coefficientsAmplitude);
            
            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                let coefficients = Mapped(ComplexPolar, Zip(coefficientsAmplitude, coefficientsPhase));
                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);
                
                // Test probability and phases of uniform superposition
                let op = StatePreparationComplexCoefficients(coefficients);
                
                using (control = Qubit[1]) {
                    
                    // Test probability
                    H(control[0]);
                    Controlled op(control, qubitsLE);
                    X(control[0]);
                    Controlled (ApplyToEachCA(H, _))(control, qubitsLE!);
                    X(control[0]);
                    
                    for ((idxCoeff, amp) in Enumerated(normalizedCoefficients)) {
                        let prob = amp * amp;
                        AssertProbInt(idxCoeff, prob, qubitsLE, tolerance);
                    }
                    
                    ResetAll(control + qubits);
                    
                    //Test phase
                    for (repeats in 0 .. nCoefficients / 2) {
                        H(control[0]);
                        Controlled op(control, qubitsLE);
                        X(control[0]);
                        Controlled (ApplyToEachCA(H, _))(control, qubitsLE!);
                        X(control[0]);
                        let indexMeasuredInteger = MeasureInteger(qubitsLE);
                        let phase = coefficientsPhase[indexMeasuredInteger];
                        Message($"StatePreparationComplexCoefficientsTest: expected phase = {phase}.");
                        Diag.AssertPhase(-0.5 * phase, control[0], tolerance);
                        ResetAll(control + qubits);
                    }
                }
            }
        }
    }
    
    
    // Test probabilities and phases of arbitrary multi-qubit superposition.
    @Diag.Test("QuantumSimulator")
    operation TestStatePreparationComplexCoefficientsArbitraryMultiQubitPhase() : Unit {
        
        let tolerance = 1E-09;

        let testCases = [
            StatePreparationTestCase(1, [1.0986553, 0.359005], [0.419893, 0.118445]),
            StatePreparationTestCase(2, [1.0986553, 0.359005, -0.123, 9.238], [0.419893, 0.118445, -0.467395, 0.419893]),
            StatePreparationTestCase(3, [1.0986553, 0.359005, 0.465689, 0.467395, 0.419893, 0.118445, 0.123, 9.238], [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445, 0.461883, 0.149609])
        ];
        
        // Loop over tests
        for (testCase in testCases) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCase!;
            let nCoefficients = Length(coefficientsAmplitude);
            
            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                let coefficients = Mapped(ComplexPolar, Zip(coefficientsAmplitude, coefficientsPhase));
                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);
                
                // Test probability and phases of arbitrary superposition
                let opComplex = StatePreparationComplexCoefficients(coefficients);
                let opReal = StatePreparationPositiveCoefficients(coefficientsAmplitude);
                
                using (control = Qubit[1]) {
                    
                    // Test probability
                    H(control[0]);
                    Controlled opComplex(control, qubitsLE);
                    X(control[0]);
                    Controlled opReal(control, qubitsLE);
                    X(control[0]);
                    
                    for ((idxCoeff, amp) in Enumerated(normalizedCoefficients)) {
                        let prob = amp * amp;
                        AssertProbInt(idxCoeff, prob, qubitsLE, tolerance);
                    }
                    
                    ResetAll(control + qubits);
                    
                    // Test phase
                    for (repeats in 0 .. nCoefficients / 2) {
                        H(control[0]);
                        Controlled opComplex(control, qubitsLE);
                        X(control[0]);
                        Controlled opReal(control, qubitsLE);
                        X(control[0]);
                        let indexMeasuredInteger = MeasureInteger(qubitsLE);
                        let phase = coefficientsPhase[indexMeasuredInteger];
                        Message($"StatePreparationComplexCoefficientsTest: expected phase = {phase}.");
                        Diag.AssertPhase(-0.5 * phase, control[0], tolerance);
                        ResetAll(control + qubits);
                    }
                }
            }
        }
    }

}


