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

    // number of qubits, abs(amplitude), phase
    newtype StatePreparationTestCase = (Int, Double[], Double[]);



    operation PrepareUniformSuperpositionTest() : Unit {
        body (...) {
            let nQubits = 5;
            using(qubits = Qubit[nQubits]) {
                for(nIndices in 1..2^nQubits)
                {
                    Message($"Testing nIndices {nIndices} on {nQubits} qubits");
                    PrepareUniformSuperposition(nIndices, LittleEndian(qubits));
                
                    ApplyToEachCA(H,qubits);

                    using(flag = Qubit[1])
                    {
                        (ControlledOnInt(0, X))(qubits, flag[0]);
                        AssertProb([PauliZ], flag, One, IntAsDouble(nIndices)/IntAsDouble(2^nQubits), "", 1e-10);
                        (ControlledOnInt(0, X))(qubits, flag[0]);
                        ApplyToEachCA(H,qubits);

                        let measuredInt = MeasureInteger(LittleEndian(qubits));

                        if(measuredInt >= nIndices){
                            fail $"Measured integer {measuredInt} which is bigger than expected of number state {nIndices}.";
                        }

                        ResetAll(flag);
                    }
                    ResetAll(qubits);
                }
            }
        }
    }

}


