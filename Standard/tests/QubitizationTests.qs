// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Testing;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    // BlockEncoding.qs tests

    // The returned operations encode the Hamiltonian (cos^2(angle) I+sin^2(angle) X)/2.
    function LCUTestHelper() : (Double[], Double, Double, (Qubit[] => Unit : Adjoint, Controlled), ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled)){
        let angle = 1.789;
        let eigenvalues = [0.5, 0.5 * Cos(angle * 2.0)];
        let prob = PowD(Cos(angle),4.0)+PowD(Sin(angle),4.0);
        let inverseAngle = ArcSin(PowD(Sin(angle),2.0)/Sqrt(prob));
        let statePreparation = Exp([PauliY], angle, _);
        let selector = Controlled (ApplyToEachCA(X, _))(_, _);
        return (eigenvalues, prob, inverseAngle, statePreparation, selector);
    }
    
    // This checks that BlockEncodingByLCU encodes the correct Hamiltonian. 
    operation BlockEncodingByLCUTest() : Unit {
        body (...) {
            let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
            let LCU = BlockEncodingByLCU(statePreparation, selector);
            using(qubits = Qubit[2]){
                let auxiliary = [qubits[0]];
                let system = [qubits[1]];

                for(rep in 0..5){
                    LCU(auxiliary, system);
                    AssertProb([PauliZ], auxiliary, Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                    let result = M(auxiliary[0]);
                    if(result == Zero) {
                        Exp([PauliY],1.0 * inverseAngle, system);
                        AssertProb([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                    }
                    ResetAll(qubits);
                }
            }
        }
    }

    // This checks that BlockEncodingReflectionByLCU encodes the correct Hamiltonian.
    operation BlockEncodingReflectionByLCUTest() : Unit {
        body (...) {
            let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
            let LCU = BlockEncodingReflectionByLCU(statePreparation, selector);
            using(qubits = Qubit[4]){
                let auxiliary = qubits[2..3];
                let system = [qubits[0]];
                let flag = qubits[1];

                for (rep in 0..5) {
                    LCU!!(auxiliary, system);
                    X(flag);
                    (ControlledOnInt(0, X))(auxiliary, flag);
                    AssertProb([PauliZ],[flag], Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                    let result = M(flag);
                    if(result == Zero) {
                        Exp([PauliY],1.0 * inverseAngle, system);
                        AssertProb([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                    }
                    ResetAll(qubits);
                }
            }
        }
    }

    // This checks that QuantumWalkByQubitization encodes the correct Hamiltonian.
    operation QuantumWalkByQubitizationTest() : Unit {
        body (...) {
            let (eigenvalues, prob, inverseAngle, statePreparation, selector) = LCUTestHelper();
            let LCU = QuantumWalkByQubitization(BlockEncodingReflectionByLCU(statePreparation, selector));
            using(qubits = Qubit[4]){
                let auxiliary = qubits[2..3];
                let system = [qubits[0]];
                let flag = qubits[1];

                for(rep in 0..5){
                    LCU(auxiliary, system);
                    X(flag);
                    (ControlledOnInt(0, X))(auxiliary, flag);
                    AssertProb([PauliZ],[flag], Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                    let result = M(flag);
                    if(result == Zero) {
                        Exp([PauliY],1.0 * inverseAngle, system);
                        AssertProb([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                    }
                    ResetAll(qubits);
                }
            }
        }
    }

    // QubitizationPauliEvolutionSet.qs tests

    // This encodes the Hamiltonian (cos^2(angle) I+sin^2(angle) X)/2.
    operation PauliBlockEncodingLCUTest() : Unit {
        body (...) {
            let angle = 0.123;
            let cosSquared = Cos(angle) * Cos(angle);
            let prob = PowD(Cos(angle),4.0)+PowD(Sin(angle),4.0);
            let inverseAngle = ArcSin(PowD(Sin(angle),2.0)/Sqrt(prob));

            mutable genIndices = new GeneratorIndex[2];
            set genIndices[0] = GeneratorIndex(([0],[cosSquared]),[0]);
            set genIndices[1] = GeneratorIndex(([1],[1.0-cosSquared]),[0]);
            
            let generatorSystem = GeneratorSystem(2, LookupFunction(genIndices));

            let (norm, LCU) = PauliBlockEncoding(generatorSystem);
            using (qubits = Qubit[2]) {
                let auxiliary = [qubits[0]];
                let system = [qubits[1]];

                for (rep in 0..5) {
                    LCU!!(auxiliary, system);
                    AssertProb([PauliZ], auxiliary, Zero, prob, "Error0: Z Success probability does not match theory", 1e-10);
                    let result = M(auxiliary[0]);
                    if(result == Zero) {
                        Exp([PauliY],1.0 * inverseAngle, system);
                        AssertProb([PauliZ], system, Zero, 1.0, "Error1: Z Success probability does not match theory", 1e-10);
                    }
                    ResetAll(qubits);
                }
            }
        }
    }

    // Array.qs tests
    function IntArrayFromRangeTest() : Unit {
        mutable testCases = new (Int[], Range)[4];
        let e = new Int[0];
        set testCases[0] = ([1, 3, 5, 7], 1..2..8);
        set testCases[1] = ([9, 6, 3, 0, -3], 9..-3..-3);
        set testCases[2] = (e, 0..2..-1);
        set testCases[3] = ([0], 0..4..3);
        for (idxTest in 0..Length(testCases) - 1) {
            let (expected, range) = testCases[idxTest];
            let output = IntArrayFromRange(range);
            Ignore(Map(AssertIntEqual(_, _, "Pad failed."), Zip(output, expected)));
        }
    }

   operation InPlaceMajorityTest() : Unit {
        body (...) {
            // Majority function truth table: x;y;z | output
            let testCases =         [[false, false, false, false],
                                     [false, false,  true, false],
                                     [false,  true, false, false],
                                     [false,  true,  true,  true],
                                     [ true, false, false, false],
                                     [ true, false,  true,  true],
                                     [ true,  true, false,  true],
                                     [ true,  true,  true,  true]];
            using (qubits = Qubit[3]) {
                for(idxTest in 0..Length(testCases)-1){
                    let testCase = testCases[idxTest];
                    Message($"Test case {idxTest}.");
                    for(idxQubit in 0..2){
                        if(testCase[idxQubit]){
                            X(qubits[idxQubit]);
                        }
                    }
                    InPlaceMajority(qubits[0], qubits[1..2]);

                    if(testCase[3] == false){
                        AssertProb([PauliZ], qubits[0..0], Zero, 1.0, "", 1e-10);
                    }
                    else{
                        AssertProb([PauliZ], qubits[0..0], One, 1.0, "", 1e-10);
                    }
                    ResetAll(qubits);
                }
            }
        }
   }

   operation ApplyRippleCarryComparatorTest() : Unit{
        body (...) {
            let nQubits = 4;
            let intMax = 2^nQubits-1;
            for(x in 0..intMax){
                for(y in 0..intMax){
                    mutable result = Zero;
                    if(x > y){
                        set result = One;
                    }
                
                    Message($"Test case. {x} > {y} = {result}");
                    using(qubits = Qubit[nQubits*2 + 1]){
                        let xRegister = BigEndian(qubits[0..nQubits-1]);
                        let yRegister = BigEndian(qubits[nQubits..2*nQubits-1]);
                        let output = qubits[2*nQubits];

                        InPlaceXorBE(x, xRegister);
                        InPlaceXorBE(y, yRegister);
                        ApplyRippleCarryComparatorBE(xRegister, yRegister, output);

                        AssertProb([PauliZ], [output], result, 1.0, "", 1e-10);
                        if(result == One){
                            X(output);
                        }

                        (Adjoint InPlaceXorBE)(y, yRegister);
                        (Adjoint InPlaceXorBE)(x, xRegister);

                    
                    }
                }
            }
        }
   }
}
