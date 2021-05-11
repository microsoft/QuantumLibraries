// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Random;

    // Tests the discretization algorithm
    @Test("QuantumSimulator")
    operation TestPurifiedMixedStateDiscretization() : Unit {
        for rep in 0..20 {
            let coeffs = DrawRandomInt(2, 5002);
            let bitsPrecision = DrawRandomInt(1, 31);
            let barHeight = 2^(bitsPrecision) - 1;
            mutable coefficients = DrawMany(DrawRandomDouble, coeffs, (0.0, 1.0));
            Message($"Test case coeffs {coeffs}, bitsPrecision {bitsPrecision}");
            // This avoids the case where coefficient are all zeros.
            let rnd = DrawRandomInt(0, coeffs - 1);
            set coefficients w/= rnd <- coefficients[rnd] + 1.0;

            let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(bitsPrecision, coefficients);

            Message($"One-norm {oneNorm}");

            // Reconstruct coefficients
            mutable coefficientsOutInt = new Int[coeffs];
            for idx in 0..coeffs - 1 {
                set coefficientsOutInt w/= idx <- coefficientsOutInt[idx] + keepCoeff[idx];
                if (altIndex[idx] >= 0) {
                    set coefficientsOutInt w/= altIndex[idx] <- coefficientsOutInt[altIndex[idx]] + barHeight - keepCoeff[idx];
                }
            }

            // Reconstruct coefficients
            mutable coefficientsOut = new Double[coeffs];
            mutable errors = new Double[coeffs];
            mutable maxError = 0.0;
            for i in 0..coeffs - 1 {
                set coefficientsOut w/= i <- oneNorm * IntAsDouble(coefficientsOutInt[i]) / IntAsDouble(barHeight * coeffs);
                let error = AbsD(coefficients[i] - coefficientsOut[i]) / oneNorm  /( PowD(2.0, IntAsDouble(-bitsPrecision)) / IntAsDouble(coeffs));
                set errors w/= i <- error;
                if (AbsD(error) > AbsD(maxError)) {
                    set maxError = error;
                }
            }
            Message($"coeffs {coeffs}, bitsPrecision {bitsPrecision}, maxError {maxError}");
            for i in 0..coeffs - 1 {
                if (errors[i] >= IntAsDouble(3)) {
                    fail $"index {i} reconstructed coefficient incorrect. Error is {errors[i]}";
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation TestPurifiedMixedState() : Unit {
        for coeffs in 2..7 {
            for nBitsPrecision in -1..-1..-2 {
                Message($"[TestPurifiedMixedState] Test case: coeffs = {coeffs}, nBitsPrecision = {nBitsPrecision}");
                let targetError = PowD(2.0, IntAsDouble(nBitsPrecision));
                let probtargetError = targetError / IntAsDouble(coeffs);
                let coefficients = DrawMany(DrawRandomDouble, coeffs, (-1.0, 1.0));

                let purifiedState = PurifiedMixedState(targetError, coefficients);

                use coeffRegister = Qubit[purifiedState::Requirements::NIndexQubits];
                use garbageQubits = Qubit[purifiedState::Requirements::NGarbageQubits];
                let coeffQubits = LittleEndian(coeffRegister);

                // Check that probability of each number state in nCoeffQubits is as expected.
                within {
                    purifiedState::Prepare(coeffQubits, new Qubit[0], garbageQubits);
                } apply {
                    for stateIndex in 0..coeffs - 1 {
                        let prob = AbsD(coefficients[stateIndex]) / purifiedState::Norm;
                        AssertProbInt(stateIndex, prob, coeffQubits, probtargetError);
                    }
                }
            }
        }
    }

    // NB: We should consider making this a public operation, perhaps part
    //     of the improvements in https://github.com/microsoft/QuantumLibraries/issues/337.
    internal operation AssertSignedProbInt(stateIndex : Int, expected : Double, sign : Qubit, qubits : LittleEndian, tolerance : Double) : Unit {
        use flag = Qubit();
        let signOffset = expected < 0.0 ? 1 <<< Length(qubits!) | 0;
        within {
            (ControlledOnInt(stateIndex + signOffset, X))(qubits! + [sign], flag);
        } apply {
            AssertMeasurementProbability([PauliZ], [flag], One, AbsD(expected), $"AssertSignedProbInt failed on stateIndex {stateIndex}, expected probability {expected}.", tolerance);
        }
    }

    @Test("QuantumSimulator")
    operation TestPurifiedMixedStateWithData() : Unit {
        for coeffs in 2..7 {
            for nBitsPrecision in -1..-1..-2 {
                let targetError = PowD(2.0, IntAsDouble(nBitsPrecision));
                let probtargetError = targetError / IntAsDouble(coeffs);
                let coefficients = DrawMany(DrawRandomDouble, coeffs, (-1.0, 1.0));
                let signs = Mapped(Compose(ConstantArray<Bool>(1, _), LessThanD(_, 0.0)), coefficients);

                let purifiedState = PurifiedMixedStateWithData(targetError, Zipped(coefficients, signs));

                use coeffRegister = Qubit[purifiedState::Requirements::NIndexQubits];
                use signQubit = Qubit();
                use garbageQubits = Qubit[purifiedState::Requirements::NGarbageQubits];
                let coeffQubits = LittleEndian(coeffRegister);

                // Check that probability of each number state in nCoeffQubits is as expected.
                within {
                    purifiedState::Prepare(coeffQubits, [signQubit], garbageQubits);
                } apply {
                    for stateIndex in 0..coeffs - 1 {
                        let prob = coefficients[stateIndex] / purifiedState::Norm;
                        AssertSignedProbInt(stateIndex, prob, signQubit, coeffQubits, probtargetError);
                    }
                }
            }
        }
    }

}
