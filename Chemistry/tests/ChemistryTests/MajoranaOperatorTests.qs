// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Testing;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    
    // Test OptimizedBEXY operator.
    operation OptimizedBEOperatorZeroTestHelper (pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int) : Unit {
        
        let indexRegisterSize = Ceiling(Lg(ToDouble(targetRegisterSize)));
        
        using (pauliBasisQubit = Qubit[1]) {
            
            using (indexRegister = Qubit[indexRegisterSize]) {
                
                using (targetRegister = Qubit[targetRegisterSize]) {
                    
                    // Choose X or Y operator.
                    if (pauliBasis == PauliX) {
                        // no op
                    }
                    elif (pauliBasis == PauliY) {
                        X(pauliBasisQubit[0]);
                    }
                    
                    // Create indexRegister state.
                    InPlaceXorBE(targetIndex, BigEndian(indexRegister));
                    
                    // Initialize targetRegister states in |0>
                    OptimizedBEXY(pauliBasisQubit[0], BigEndian(indexRegister), targetRegister);
                    
                    for (idxTest in 0 .. targetRegisterSize - 1) {
                        let testQubit = targetRegister[idxTest];
                        
                        if (targetIndex > idxTest) {
                            
                            // Test Z Pauli
                            // |0> -> |0>
                            // |+> -> |->
                            Message($"Test Z Pauli on qubit {idxTest}");
                            AssertProb([PauliZ], [testQubit], Zero, 1.0, $"Error: Test {idxTest} {idxTest} Z Pauli |0>", 1E-10);
                        }
                        elif (targetIndex == idxTest) {
                            
                            // Test X Pauli
                            // |0> -> |1>
                            // |+> -> |+>
                            
                            // Test Y Pauli
                            // |0> -> i|1>
                            // |+> -> -i|->
                            Message($"Test X or Y Pauli on qubit {idxTest}");
                            AssertProb([PauliZ], [testQubit], One, 1.0, $"Error: Test {idxTest} X or Y Pauli |0>", 1E-10);
                        }
                        else {
                            
                            // Test Identitfy Pauli
                            // |0> -> |0>
                            // |+> -> |+>
                            Message($"Test ZI Pauli on qubit {idxTest}");
                            AssertProb([PauliZ], [testQubit], Zero, 1.0, $"Error: Test {idxTest} I Pauli |0>", 1E-10);
                        }
                    }
                    
                    OptimizedBEXY(pauliBasisQubit[0], BigEndian(indexRegister), targetRegister);
                    Adjoint InPlaceXorBE(targetIndex, BigEndian(indexRegister));
                    
                    // Choose X or Y operator.
                    if (pauliBasis == PauliX) {
                        // no op
                    }
                    elif (pauliBasis == PauliY) {
                        X(pauliBasisQubit[0]);
                    }
                }
            }
        }
    }
    
    
    // Test OptimizedBEXY operator.
    operation OptimizedBEOperatorPlusTestHelper (pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int) : Unit {
        
        let indexRegisterSize = Ceiling(Lg(ToDouble(targetRegisterSize)));
        using(pauliBasisQubit = Qubit[1]){
            using(indexRegister = Qubit[indexRegisterSize]){
                using(targetRegister = Qubit[targetRegisterSize]){            
                    // Choose X or Y operator.
                    if(pauliBasis == PauliX){
                        // no op
                    }
                    elif(pauliBasis == PauliY){
                        X(pauliBasisQubit[0]);
                    }

                    // Create indexRegister state.
                    InPlaceXorBE(targetIndex, BigEndian(indexRegister));
                        
                    // Initialize targetRegister states in |+>
                    ApplyToEachCA(H, targetRegister);
                    OptimizedBEXY(pauliBasisQubit[0], BigEndian(indexRegister), targetRegister);
                    for(idxTest in 0..targetRegisterSize-1){
                        let testQubit = targetRegister[idxTest];
                        if(targetIndex > idxTest){
                            // Test Z Pauli 
                            // |0> -> |0>
                            // |+> -> |->
                            Message($"Test Z Pauli on qubit {idxTest}");
                            AssertProb([PauliX], [testQubit], One, 1.0, $"Error: Test {idxTest} Z Pauli |+>", 1e-10);
                        }
                        elif(targetIndex == idxTest){
                            // Test X Pauli
                            // |0> -> |1>
                            // |+> -> |+>
                            if(pauliBasis == PauliX){
                                Message($"Test X Pauli on qubit {idxTest}");
                                AssertProb([PauliX], [testQubit], Zero, 1.0, $"Error: Test {idxTest} X Pauli |+>", 1e-10);
                            }
                                
                            // Test Y Pauli
                            // |0> -> i|1>
                            // |+> -> -i|->
                            if(pauliBasis == PauliY){
                                Message($"Test Y Pauli on qubit {idxTest}");
                                AssertProb([PauliX], [testQubit], One, 1.0, $"Error: Test {idxTest} Y Pauli |+>", 1e-10);
                            }

                        }
                        else{
                            // Test Identitfy Pauli
                            // |0> -> |0>
                            // |+> -> |+>
                            Message($"Test I Pauli on qubit {idxTest}");
                            AssertProb([PauliX], [testQubit], Zero, 1.0, $"Error: Test {idxTest} I Pauli |+>", 1e-10);
                        }
                    }
                    OptimizedBEXY(pauliBasisQubit[0], BigEndian(indexRegister), targetRegister);
                    ApplyToEachCA(H, targetRegister);

                    (Adjoint InPlaceXorBE)(targetIndex, BigEndian(indexRegister));

                    // Choose X or Y operator.
                    if(pauliBasis == PauliX){
                        // no op
                    }
                    elif(pauliBasis == PauliY){
                        X(pauliBasisQubit[0]);
                    }
                }
            }
        }
    }
    
    
    operation OptimizedBEOperatorZeroTest () : Unit {
        
        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;
        
        for (idxPauli in 0 .. 1) {
            let pauliBasis = paulis[idxPauli];
            
            for (targetIndex in 0 .. targetRegisterSize - 1) {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                OptimizedBEOperatorZeroTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }
    
    
    operation OptimizedBEOperatorPlusTest () : Unit {
        
        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;
        
        for (idxPauli in 0 .. 1) {
            let pauliBasis = paulis[idxPauli];
            
            for (targetIndex in 0 .. targetRegisterSize - 1) {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                OptimizedBEOperatorPlusTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }
    
    
    // Test phase of controlled OptimizedBEXY operator.
    operation ControlledOptimizedBEOperatorTestHelper (pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int) : Unit {
        
        let indexRegisterSize = Ceiling(Lg(ToDouble(targetRegisterSize)));
        
        using (pauliBasisQubit = Qubit[1]) {
            
            using (indexRegister = Qubit[indexRegisterSize]) {
                
                using (targetRegister = Qubit[targetRegisterSize]) {
                    
                    using (controlRegister = Qubit[1]) {
                        let testQubit = targetRegister[targetIndex];
                        
                        // Create indexRegister state.
                        InPlaceXorBE(targetIndex, BigEndian(indexRegister));
                        
                        // Initialize control in |+> state.
                        H(controlRegister[0]);
                        
                        // Choose X or Y operator.
                        if (pauliBasis == PauliX) {
                            
                            // Initialize testQubit state in X +1 eigenstate
                            H(testQubit);
                            Controlled OptimizedBEXY(controlRegister, (pauliBasisQubit[0], BigEndian(indexRegister), targetRegister));
                            AssertPhase(0.0, controlRegister[0], 1E-10);
                            Adjoint Controlled OptimizedBEXY(controlRegister, (pauliBasisQubit[0], BigEndian(indexRegister), targetRegister));
                            H(testQubit);
                        }
                        elif (pauliBasis == PauliY) {
                            X(pauliBasisQubit[0]);
                            
                            // Initialize testQubit state Y +1 eigenstate
                            H(testQubit);
                            S(testQubit);
                            Controlled OptimizedBEXY(controlRegister, (pauliBasisQubit[0], BigEndian(indexRegister), targetRegister));
                            AssertPhase(0.0, controlRegister[0], 1E-10);
                            Adjoint Controlled OptimizedBEXY(controlRegister, (pauliBasisQubit[0], BigEndian(indexRegister), targetRegister));
                            Adjoint S(testQubit);
                            H(testQubit);
                            X(pauliBasisQubit[0]);
                        }
                        
                        H(controlRegister[0]);
                        Adjoint InPlaceXorBE(targetIndex, BigEndian(indexRegister));
                    }
                }
            }
        }
    }
    
    
    operation ControlledOptimizedBEOperatorPlusTest () : Unit {
        
        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;
        
        for (idxPauli in 0 .. 1) {
            let pauliBasis = paulis[idxPauli];
            
            for (targetIndex in 0 .. targetRegisterSize - 1) {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                ControlledOptimizedBEOperatorTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }
    
    
    // Test SelectZ operator
    operation SelectZTest () : Unit {
        
        let targetRegisterSize = 7;
        let indexRegisterSize = Ceiling(Lg(ToDouble(targetRegisterSize)));
        
        using (targetRegister = Qubit[targetRegisterSize]) {
            
            using (indexRegister = Qubit[indexRegisterSize]) {
                
                for (idxTest in 0 .. targetRegisterSize - 1) {
                    H(targetRegister[idxTest]);
                    InPlaceXorLE(idxTest, LittleEndian(Reverse(indexRegister)));
                    SelectZ(BigEndian(indexRegister), targetRegister);
                    AssertProb([PauliX], [targetRegister[idxTest]], One, 1.0, $"Error: Test {idxTest} X Pauli |+>", 1E-10);
                    Z(targetRegister[idxTest]);
                    Adjoint InPlaceXorLE(idxTest, LittleEndian(Reverse(indexRegister)));
                    H(targetRegister[idxTest]);
                }
            }
        }
    }
    
}


