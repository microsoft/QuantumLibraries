// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    // Tests the discretization algorithm
    operation QuantumROMDiscretizationTest() : Unit {
        for(rep in 0..20){
            let coeffs = RandomInt(5000)+2;
            let bitsPrecision = RandomInt(30)+1;
            let barHeight = 2^(bitsPrecision) - 1;
            mutable coefficients = new Double[coeffs];
            Message($"Test case coeffs {coeffs}, bitsPrecision {bitsPrecision}");
            for (idx in 0..coeffs-1)
            {
                set coefficients w/= idx <- 1000.0 * RandomReal(2*bitsPrecision);
            }
            // This avoids the case where coefficient are all zeros.
            let rnd = RandomInt(coeffs);
            set coefficients w/= rnd <- coefficients[rnd] + 1.0;

            let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(bitsPrecision, coefficients);

            Message($"One-norm {oneNorm}");

            // Reconstruct coefficients
            mutable coefficientsOutInt = new Int[coeffs];
            for (idx in 0..coeffs-1)
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
            for (i in 0..coeffs-1)
            {
                set coefficientsOut w/= i <- oneNorm * IntAsDouble(coefficientsOutInt[i]) / IntAsDouble(barHeight * coeffs);
                let error = AbsD(coefficients[i] - coefficientsOut[i]) / oneNorm  /( PowD(2.0, IntAsDouble(-bitsPrecision)) / IntAsDouble(coeffs));
                set errors w/= i <- error;
                if(AbsD(error) > AbsD(maxError)){
                    set maxError = error;
                }
            }
            Message($"coeffs {coeffs}, bitsPrecision {bitsPrecision}, maxError {maxError}");
            for(i in 0..coeffs-1){
                if(errors[i] < IntAsDouble(3)){
                    // test passes
                }
                else{
                    fail $"index {i} reconstructed coefficient incorrect. Error is {errors[i]}";
                }
            }
        }
    }

    operation QuantumROMTest() : Unit {
        for(coeffs in 2..7){
            for(nBitsPrecision in -1..-1..-2){
                let targetError = PowD(2.0, IntAsDouble(nBitsPrecision));
                let probtargetError = targetError / IntAsDouble(coeffs);
                mutable coefficients = new Double[coeffs];
                for (idx in 0..coeffs-1)
                {
                    set coefficients w/= idx <- RandomReal(2*32);
                }
                let ((nTotal, (nCoeffQubits, nGarbageQubits)), oneNorm, op) =  QuantumROM(targetError, coefficients);
                Message($"Test case coeffs {coeffs}, bitsPrecision {nCoeffQubits}, global targetError {targetError}, probability error {probtargetError}.");
                for (idx in 0..coeffs-1)
                {
                    let tmp = AbsD(coefficients[idx]) / oneNorm;
                    Message($"{idx} expected prob = {tmp}.");
                }

                Message($"Qubits used: {nGarbageQubits} + {nCoeffQubits}");
                using(qubits = Qubit[nTotal]){
                    let (register, rest) = _QuantumROMQubitManager(targetError, coeffs, qubits);
                    let (coeffQubits, garbageQubits) = register;
                    op(register);
                    // Now check that probability of each number state in nCoeffQubits is as expected.
                    for(stateIndex in 0..coeffs-1){
                        let prob = AbsD(coefficients[stateIndex]) / oneNorm;
                        Message($"Testing probability {prob} on index {stateIndex}");
                        //BAssertProbIntBE(stateIndex, AbsD(coefficients[stateIndex]) / oneNorm, BigEndian(coeffQubits), targetError / IntAsDouble(coeffs));
                    }

                    (Adjoint op)(register);

                }
            }
        }
    }

}
