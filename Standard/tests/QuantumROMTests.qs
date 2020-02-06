// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    // Tests the discretization algorithm
    @Test("QuantumSimulator")
    operation _QuantumROMDiscretizationTest() : Unit {
        for (rep in 0..20) {
            let nCoefficients = RandomInt(5000) + 2;
            let nBitsPrecision = RandomInt(30) + 1;
            let barHeight = 2 ^ nBitsPrecision - 1;
            mutable coefficients = new Double[nCoefficients];
            Message($"Test case coeffs {nCoefficients}, bitsPrecision {nBitsPrecision}");
            for (idx in 0..nCoefficients - 1) {
                set coefficients w/= idx <- 1000.0 * RandomReal(2 * nBitsPrecision);
            }
            // This avoids the case where coefficient are all zeros.
            let rnd = RandomInt(nCoefficients);
            set coefficients w/= rnd <- coefficients[rnd] + 1.0;

            let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, coefficients);

            Message($"One-norm {oneNorm}");

            // Reconstruct coefficients
            mutable coefficientsOutInt = new Int[nCoefficients];
            for (idx in 0..nCoefficients - 1) {
                set coefficientsOutInt w/= idx <- coefficientsOutInt[idx] + keepCoeff[idx];
                if (altIndex[idx] >= 0) {
                    set coefficientsOutInt w/= altIndex[idx] <- coefficientsOutInt[altIndex[idx]] + barHeight - keepCoeff[idx];
                }
            }

            // Reconstruct coefficients
            mutable coefficientsOut = new Double[nCoefficients];
            mutable errors = new Double[nCoefficients];
            mutable maxError = 0.0;
            for (i in 0..nCoefficients - 1) {
                set coefficientsOut w/= i <- oneNorm * IntAsDouble(coefficientsOutInt[i]) / IntAsDouble(barHeight * nCoefficients);
                let error = AbsD(coefficients[i] - coefficientsOut[i]) / oneNorm  /( PowD(2.0, IntAsDouble(-nBitsPrecision)) / IntAsDouble(nCoefficients));
                set errors w/= i <- error;
                set maxError = MaxD(AbsD(error), maxError);
            }
            Message($"coeffs {nCoefficients}, bitsPrecision {nBitsPrecision}, maxError {maxError}");
            for (i in 0..nCoefficients - 1) {
                if (errors[i] >= 3.0) {
                    fail $"index {i} reconstructed coefficient incorrect. Error is {errors[i]}";
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation QuantumROMTest() : Unit {
        for (coeffs in 2..7) {
            for (nBitsPrecision in -1..-1..-2) {
                let targetError = PowD(2.0, IntAsDouble(nBitsPrecision));
                let probtargetError = targetError / IntAsDouble(coeffs);
                let coefficients = ForEach(Delay(RandomReal, 2 * 32, _), ConstantArray(coeffs, ()));
                let preparation = PurifiedMixedState(targetError, coefficients);
                Message($"Test case coeffs {coeffs}, bitsPrecision {preparation::Requirements::NIndexQubits}, global targetError {targetError}, probability error {probtargetError}.");
                for (idx in 0..coeffs - 1) {
                    let tmp = AbsD(coefficients[idx]) / preparation::Norm;
                    Message($"{idx} expected prob = {tmp}.");
                }

                Message($"Qubits used: {preparation::Requirements::NGarbageQubits} + {preparation::Requirements::NIndexQubits}");
                using (qubits = Qubit[preparation::Requirements::NTotalQubits]) {
                    let (register, rest) = _PartitionedForQuantumROM(targetError, coeffs, qubits);
                    let (coeffQubits, garbageQubits) = register;
                    within {
                        preparation::Prepare(register);
                    } apply {
                        // Now check that probability of each number state in nCoeffQubits is as expected.
                        for (stateIndex in 0..coeffs - 1) {
                            let prob = AbsD(coefficients[stateIndex]) / preparation::Norm;
                            Message($"Testing probability {prob} on index {stateIndex}");
                            AssertProbInt(
                                stateIndex,
                                prob,
                                coeffQubits,
                                targetError / IntAsDouble(coeffs)
                            );
                        }
                    }

                }
            }
        }
    }

}
