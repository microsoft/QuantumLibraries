// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Tests {
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
    operation TestQuantumROMDiscretization() : Unit {
        for(rep in 0..20){
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
            for (idx in 0..coeffs - 1)
            {
                set coefficientsOutInt w/= idx <- coefficientsOutInt[idx] + keepCoeff[idx];
                if (altIndex[idx] >= 0)
                {
                    set coefficientsOutInt w/= altIndex[idx] <- coefficientsOutInt[altIndex[idx]] + barHeight - keepCoeff[idx];
                }
            }

            // Reconstruct coefficients
            mutable coefficientsOut = new Double[coeffs];
            mutable errors = new Double[coeffs];
            mutable maxError = 0.0;
            for (i in 0..coeffs - 1) {
                set coefficientsOut w/= i <- oneNorm * IntAsDouble(coefficientsOutInt[i]) / IntAsDouble(barHeight * coeffs);
                let error = AbsD(coefficients[i] - coefficientsOut[i]) / oneNorm  /( PowD(2.0, IntAsDouble(-bitsPrecision)) / IntAsDouble(coeffs));
                set errors w/= i <- error;
                if (AbsD(error) > AbsD(maxError)) {
                    set maxError = error;
                }
            }
            Message($"coeffs {coeffs}, bitsPrecision {bitsPrecision}, maxError {maxError}");
            for (i in 0..coeffs - 1) {
                if (errors[i] >= IntAsDouble(3)) {
                    fail $"index {i} reconstructed coefficient incorrect. Error is {errors[i]}";
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation TestQuantumROM() : Unit {
        for(coeffs in 2..7){
            for(nBitsPrecision in -1..-1..-2){
                let targetError = PowD(2.0, IntAsDouble(nBitsPrecision));
                let probtargetError = targetError / IntAsDouble(coeffs);
                let coefficients = DrawMany(DrawRandomDouble, coeffs, (-1.0, 1.0));

                if (true) { // quantum ROM without sign
                    let ((nTotal, (nCoeffQubits, nGarbageQubits)), oneNorm, op) =  QuantumROM(targetError, coefficients);

                    using ((coeffRegister, garbageQubits) = (Qubit[nCoeffQubits], Qubit[nGarbageQubits])) {
                        let coeffQubits = LittleEndian(coeffRegister);

                        // Check that probability of each number state in nCoeffQubits is as expected.
                        within {
                            op(coeffQubits, garbageQubits);
                        } apply {
                            for (stateIndex in 0..coeffs - 1) {
                                let prob = AbsD(coefficients[stateIndex]) / oneNorm;
                                AssertProbInt(stateIndex, prob, coeffQubits, probtargetError);
                            }
                        }
                    }
                }

                if (true) { // quantum ROM with sign
                    let ((nTotal, (nCoeffQubits, nGarbageQubits)), oneNorm, op) =  QuantumROMWithSign(targetError, coefficients);

                    using ((coeffRegister, signQubit, garbageQubits) = (Qubit[nCoeffQubits], Qubit(), Qubit[nGarbageQubits])) {
                        let coeffQubits = LittleEndian(coeffRegister);

                        // Check that probability of each number state in nCoeffQubits is as expected.
                        within {
                            op(coeffQubits, signQubit, garbageQubits);
                        } apply {
                            for (stateIndex in 0..coeffs - 1) {
                                let prob = AbsD(coefficients[stateIndex]) / oneNorm;
                                AssertProbInt(stateIndex, prob, coeffQubits, probtargetError);
                            }
                        }
                    }
                }
            }
        }
    }

}
