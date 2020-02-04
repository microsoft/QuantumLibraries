// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    // number of qubits, abs(amplitude), phase
    newtype StatePreparationTestCase = (
        NQubits: Int,
        Magnitudes: Double[],
        Phases: Double[]
    );

    @Test("QuantumSimulator")
    operation StatePreparationPositiveCoefficientsTest () : Unit {
        let tolerance = 1E-09;
        let testCases = [
            // Test positive coefficients.
            StatePreparationTestCase(1, [0.773761, 0.633478], [0.0, 0.0]),
            StatePreparationTestCase(2, [0.183017, 0.406973, 0.604925, 0.659502], [0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(3, [0.0986553, 0.359005, 0.465689, 0.467395, 0.419893, 0.118445, 0.461883, 0.149609], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            StatePreparationTestCase(4, [0.271471, 0.0583654, 0.11639, 0.36112, 0.307383, 0.193371, 0.274151, 0.332542, 0.130172, 0.222546, 0.314879, 0.210704, 0.212429, 0.245518, 0.30666, 0.22773], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),

            // Test negative coefficients. Should give same probabilities as positive coefficients.
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
        for ((idxTestCase, testCase) in Enumerated(testCases)) {

            // Test negative coefficients. Should give same results as positive coefficients.
            using (qubits = Qubit[testCase::NQubits]) {
                let qubitsLE = LittleEndian(qubits);
                let op = StatePreparationPositiveCoefficients(testCase::Magnitudes);
                op(qubitsLE);
                let normalizedCoefficients = PNormalized(2.0, testCase::Magnitudes);

                for ((idxCoefficient, coefficient) in Enumerated(normalizedCoefficients)) {
                    AssertProbInt(idxCoefficient, coefficient * coefficient, qubitsLE, tolerance);
                }

                ResetAll(qubits);
            }
        }
    }


    // Test phase factor on 1-qubit uniform superposition.
    @Test("QuantumSimulator")
    operation StatePreparationComplexCoefficientsQubitPhaseTest () : Unit {

        let tolerance = 1E-09;
        mutable testCases = new StatePreparationTestCase[10];
        mutable nTests = 0;

        // Test phase factor on uniform superposition.
        set testCases w/= nTests <- StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.01]);
        set nTests = nTests + 1;
        set testCases w/= nTests <- StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.05]);
        set nTests = nTests + 1;

        // Loop over tests
        for (idxTest in 0 .. nTests - 1) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCases[idxTest]!;
            Message($"Test case {idxTest}");
            let nCoefficients = Length(coefficientsAmplitude);

            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                mutable coefficients = new ComplexPolar[nCoefficients];
                mutable coefficientsPositive = new Double[nCoefficients];

                for (idxCoeff in 0 .. nCoefficients - 1) {
                    set coefficients w/= idxCoeff <- ComplexPolar(coefficientsAmplitude[idxCoeff], coefficientsPhase[idxCoeff]);
                    set coefficientsPositive w/= idxCoeff <- coefficientsAmplitude[idxCoeff];
                }

                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);

                // Test phase factor on uniform superposition
                let phase = 0.5 * (coefficientsPhase[0] - coefficientsPhase[1]);
                let amp = normalizedCoefficients[0];
                let prob = amp * amp;
                let op = StatePreparationComplexCoefficients(coefficients);
                op(qubitsLE);
                AssertProbInt(0, prob, qubitsLE, tolerance);
                AssertProbInt(1, prob, qubitsLE, tolerance);
                AssertPhase(phase, (qubitsLE!)[0], tolerance);
                ResetAll(qubits);
            }
        }
    }


    // Test probabilities and phases factor of multi-qubit uniform superposition.
    @Test("QuantumSimulator")
    operation StatePreparationComplexCoefficientsMultiQubitPhaseTest () : Unit {

        let tolerance = 1E-09;
        mutable testCases = new StatePreparationTestCase[10];
        mutable nTests = 0;

        // Test probability and phases of uniform superposition.
        set testCases w/= nTests <- StatePreparationTestCase(1, [1.0, 1.0], [0.01, -0.01]);
        set nTests = nTests + 1;
        set testCases w/= nTests <- StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        set nTests = nTests + 1;
        set testCases w/= nTests <- StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], ConstantArray(8, PI()));
        set nTests = nTests + 1;
        set testCases w/= nTests <- StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01]);
        set nTests = nTests + 1;
        set testCases w/= nTests <- StatePreparationTestCase(3, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445, 0.461883, 0.149609]);
        set nTests = nTests + 1;

        // Loop over tests
        for (idxTest in 0 .. nTests - 1) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCases[idxTest]!;
            Message($"Test case {idxTest}");
            let nCoefficients = Length(coefficientsAmplitude);

            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                mutable coefficients = new ComplexPolar[nCoefficients];
                mutable coefficientsPositive = new Double[nCoefficients];

                for (idxCoeff in 0 .. nCoefficients - 1) {
                    set coefficients w/= idxCoeff <- ComplexPolar(coefficientsAmplitude[idxCoeff], coefficientsPhase[idxCoeff]);
                    set coefficientsPositive w/= idxCoeff <- coefficientsAmplitude[idxCoeff];
                }

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

                    for (idxCoeff in 0 .. nCoefficients - 1) {
                        let amp = normalizedCoefficients[idxCoeff];
                        let prob = amp * amp;
                        AssertProbInt(idxCoeff, prob, qubitsLE, tolerance);
                    }

                    ResetAll(control);
                    ResetAll(qubits);

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
                        AssertPhase(-0.5 * phase, control[0], tolerance);
                        ResetAll(control);
                        ResetAll(qubits);
                    }
                }
            }
        }
    }


    // Test probabilities and phases of arbitrary multi-qubit superposition.
    @Test("QuantumSimulator")
    operation StatePreparationComplexCoefficientsArbitraryMultiQubitPhaseTest () : Unit {

        let tolerance = 1E-09;
        let testCases = [
            StatePreparationTestCase(1, [1.0986553, 0.359005], [0.419893, 0.118445]),
            StatePreparationTestCase(2, [1.0986553, 0.359005, -0.123, 9.238], [0.419893, 0.118445, -0.467395, 0.419893]),StatePreparationTestCase(3, [1.0986553, 0.359005, 0.465689, 0.467395, 0.419893, 0.118445, 0.123, 9.238], [1.0986553, 0.359005, 0.465689, -0.467395, 0.419893, 0.118445, 0.461883, 0.149609])
        ];

        // Loop over tests
        for ((idxTestCase, testCase) in Enumerated(testCases)) {
            let (nQubits, coefficientsAmplitude, coefficientsPhase) = testCase!;
            Message($"Test case {idxTestCase}");
            let nCoefficients = Length(coefficientsAmplitude);

            using (qubits = Qubit[nQubits]) {
                let qubitsLE = LittleEndian(qubits);
                mutable coefficientsPositive = new Double[nCoefficients];

                let coefficients = Mapped(ComplexPolar, Zip(coefficientsAmplitude, coefficientsPhase));
                let normalizedCoefficients = PNormalized(2.0, coefficientsAmplitude);

                // Test probability and phases of arbitrary superposition
                let opComplex = StatePreparationComplexCoefficients(coefficients);
                let opReal = StatePreparationPositiveCoefficients(coefficientsAmplitude);

                using (control = Qubit()) {

                    // Test probability
                    H(control);
                    Controlled opComplex([control], qubitsLE);
                    within {
                        X(control);
                    } apply {
                        Controlled opReal([control], qubitsLE);
                    }

                    for ((idxCoeff, coeff) in Enumerated(normalizedCoefficients)) {
                        AssertProbInt(idxCoeff, coeff * coeff, qubitsLE, tolerance);
                    }

                    Reset(control);
                    ResetAll(qubits);

                    // Test phase
                    for (repeats in 0 .. nCoefficients / 2) {
                        H(control);
                        Controlled opComplex([control], qubitsLE);
                        within {
                            X(control);
                        } apply {
                            Controlled opReal([control], qubitsLE);
                        }
                        let indexMeasuredInteger = MeasureInteger(qubitsLE);
                        let phase = coefficientsPhase[indexMeasuredInteger];
                        Message($"StatePreparationComplexCoefficientsTest: expected phase = {phase}.");
                        AssertPhase(-0.5 * phase, control, tolerance);
                        Reset(control);
                        ResetAll(qubits);
                    }
                }
            }
        }
    }

}


