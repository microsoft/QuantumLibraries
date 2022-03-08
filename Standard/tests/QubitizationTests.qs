// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    // BlockEncoding.qs tests

    // The returned operations encode the Hamiltonian (cos^2(angle) I+sin^2(angle) X)/2.
    function LCUTestHelper() : (Double[], Double, Double, (Qubit[] => Unit is Adj + Ctl), ((Qubit[], Qubit[]) => Unit is Adj + Ctl)){
        let angle = 1.789;
        let eigenvalues = [0.5, 0.5 * Cos(angle * 2.0)];
        let prob = PowD(Cos(angle),4.0)+PowD(Sin(angle),4.0);
        let inverseAngle = ArcSin(PowD(Sin(angle),2.0)/Sqrt(prob));
        let statePreparation = Exp([PauliY], angle, _);
        let selector = Controlled (ApplyToEachCA(X, _))(_, _);
        return (eigenvalues, prob, inverseAngle, statePreparation, selector);
    }

    // This checks that BlockEncodingByLCU encodes the correct Hamiltonian.
    @Test("QuantumSimulator")
    operation TestBlockEncodingByLCU() : Unit {
        let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
        let LCU = BlockEncodingByLCU(statePreparation, selector);
        use qubits = Qubit[2];
        let auxiliary = [qubits[0]];
        let system = [qubits[1]];

        for rep in 0..5 {
            LCU(auxiliary, system);
            AssertMeasurementProbability([PauliZ], auxiliary, Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
            let result = M(auxiliary[0]);
            if(result == Zero) {
                Exp([PauliY],1.0 * inverseAngle, system);
                AssertMeasurementProbability([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
            }
            ResetAll(qubits);
        }
    }

    // This checks that BlockEncodingReflectionByLCU encodes the correct Hamiltonian.
    @Test("QuantumSimulator")
    operation TestBlockEncodingReflectionByLCU() : Unit {
        body (...) {
            let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
            let LCU = BlockEncodingReflectionByLCU(statePreparation, selector);
            use flag = Qubit();
            use system = Qubit[1];
            use auxiliary = Qubit[2];

            for rep in 0..5 {
                LCU!!(auxiliary, system);
                X(flag);
                (ControlledOnInt(0, X))(auxiliary, flag);
                AssertMeasurementProbability([PauliZ],[flag], Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                let result = M(flag);
                if result == Zero {
                    Exp([PauliY], 1.0 * inverseAngle, system);
                    AssertMeasurementProbability([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                }
                ResetAll([flag] + system + auxiliary);
            }
        }
    }

    // This checks that QuantumWalkByQubitization encodes the correct Hamiltonian.
    @Test("QuantumSimulator")
    operation TestQuantumWalkByQubitization() : Unit {
        body (...) {
            let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
            let LCU = QuantumWalkByQubitization(BlockEncodingReflectionByLCU(statePreparation, selector));
            use flag = Qubit();
            use system = Qubit[1];
            use auxiliary = Qubit[2];

            for rep in 0..5 {
                LCU(auxiliary, system);
                X(flag);
                (ControlledOnInt(0, X))(auxiliary, flag);
                AssertMeasurementProbability([PauliZ],[flag], Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                let result = M(flag);
                if(result == Zero) {
                    Exp([PauliY],1.0 * inverseAngle, system);
                    AssertMeasurementProbability([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                }
                ResetAll(system);
                Reset(flag);
                ResetAll(auxiliary);
            }
        }
    }

    // QubitizationPauliEvolutionSet.qs tests

    // This encodes the Hamiltonian (cos^2(angle) I+sin^2(angle) X)/2.
    @Test("QuantumSimulator")
    operation TestPauliBlockEncodingLCU() : Unit {
        body (...) {
            let angle = 0.123;
            let cosSquared = Cos(angle) * Cos(angle);
            let prob = PowD(Cos(angle),4.0)+PowD(Sin(angle),4.0);
            let inverseAngle = ArcSin(PowD(Sin(angle),2.0)/Sqrt(prob));

            let genIndices = [
                GeneratorIndex(([0],[cosSquared]),[0]),
                GeneratorIndex(([1],[1.0-cosSquared]),[0])
            ];
            
            let generatorSystem = GeneratorSystem(2, LookupFunction(genIndices));

            let (norm, LCU) = PauliBlockEncoding(generatorSystem);
            use qubits = Qubit[2];
            let auxiliary = [qubits[0]];
            let system = [qubits[1]];

            for rep in 0..5 {
                LCU!!(auxiliary, system);
                AssertMeasurementProbability([PauliZ], auxiliary, Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                let result = M(auxiliary[0]);
                if(result == Zero) {
                    Exp([PauliY],1.0 * inverseAngle, system);
                    AssertMeasurementProbability([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                }
                ResetAll(qubits);
            }
        }
    }

    // Array.qs tests
    @Test("QuantumSimulator")
    function TestRangeAsIntArray() : Unit {
        let testCases = [
            ([1, 3, 5, 7], 1..2..8),
            ([9, 6, 3, 0, -3], 9..-3..-3),
            ([], 0..2..-1),
            ([0], 0..4..3)
        ];
        for (expected, range) in testCases {
            let output = RangeAsIntArray(range);
            AllEqualityFactI(output, expected, "RangeAsIntArray failed");
        }
    }

   @Test("QuantumSimulator")
   operation TestInPlaceMajority() : Unit {
        // Majority function truth table: x;y;z | output
        let testCases =         [[false, false, false, false],
                                    [false, false,  true, false],
                                    [false,  true, false, false],
                                    [false,  true,  true,  true],
                                    [ true, false, false, false],
                                    [ true, false,  true,  true],
                                    [ true,  true, false,  true],
                                    [ true,  true,  true,  true]];
        use qubits = Qubit[3];
        for idxTest in IndexRange(testCases) {
            let testCase = testCases[idxTest];
            Message($"Test case {idxTest}.");
            for idxQubit in 0..2 {
                if(testCase[idxQubit]){
                    X(qubits[idxQubit]);
                }
            }
            ApplyMajorityInPlace(qubits[0], qubits[1..2]);

            if (testCase[3] == false) {
                AssertMeasurementProbability([PauliZ], qubits[0..0], Zero, 1.0, "", 1e-10);
            }
            else {
                AssertMeasurementProbability([PauliZ], qubits[0..0], One, 1.0, "", 1e-10);
            }
            ResetAll(qubits);
        }
   }

    @Test("QuantumSimulator")
    operation TestApplyRippleCarryComparator() : Unit {
        let nQubits = 4;
        let intMax = 2 ^ nQubits - 1;
        for x in 0..intMax {
            for y in 0..intMax {
                let result = x > y ? One | Zero;

                Message($"Test case. {x} > {y} = {result}");
                use qubits = Qubit[nQubits * 2 + 1];
                let xRegister = LittleEndian(qubits[0..nQubits-1]);
                let yRegister = LittleEndian(qubits[nQubits..2*nQubits-1]);
                let output = qubits[2*nQubits];

                ApplyXorInPlace(x, xRegister);
                ApplyXorInPlace(y, yRegister);
                CompareUsingRippleCarry(xRegister, yRegister, output);

                AssertMeasurementProbability([PauliZ], [output], result, 1.0, "", 1e-10);
                if result == One {
                    X(output);
                }

                (Adjoint ApplyXorInPlace)(y, yRegister);
                (Adjoint ApplyXorInPlace)(x, xRegister);
            }
        }
    }
}
