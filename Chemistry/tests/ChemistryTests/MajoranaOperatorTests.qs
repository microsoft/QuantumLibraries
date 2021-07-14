// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics as Diag;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

    // Test OptimizedBEXY operator.
    operation OptimizedBEOperatorZeroTestHelper (pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int) : Unit {
        let indexRegisterSize = Ceiling(Lg(IntAsDouble(targetRegisterSize)));
        use pauliBasisQubit = Qubit();
        use indexRegister = Qubit[indexRegisterSize];
        use targetRegister = Qubit[targetRegisterSize];

        // Choose X or Y operator.
        if pauliBasis == PauliX {
            // no op
        } elif pauliBasis == PauliY {
            X(pauliBasisQubit);
        }

        // Create indexRegister state.
        ApplyXorInPlace(targetIndex, LittleEndian(indexRegister));

        // Initialize targetRegister states in |0>
        OptimizedBEXY(pauliBasisQubit, LittleEndian(indexRegister), targetRegister);

        for idxTest in 0 .. targetRegisterSize - 1 {
            let testQubit = targetRegister[idxTest];

            if targetIndex > idxTest {

                // Test Z Pauli
                // |0> -> |0>
                // |+> -> |->
                Message($"Test Z Pauli on qubit {idxTest}");
                Diag.AssertMeasurementProbability([PauliZ], [testQubit], Zero, 1.0, $"Error: Test {idxTest} {idxTest} Z Pauli |0>", 1E-10);
            } elif targetIndex == idxTest {

                // Test X Pauli
                // |0> -> |1>
                // |+> -> |+>

                // Test Y Pauli
                // |0> -> i|1>
                // |+> -> -i|->
                Message($"Test X or Y Pauli on qubit {idxTest}");
                Diag.AssertMeasurementProbability([PauliZ], [testQubit], One, 1.0, $"Error: Test {idxTest} X or Y Pauli |0>", 1E-10);
            } else {

                // Test Identitfy Pauli
                // |0> -> |0>
                // |+> -> |+>
                Message($"Test ZI Pauli on qubit {idxTest}");
                Diag.AssertMeasurementProbability([PauliZ], [testQubit], Zero, 1.0, $"Error: Test {idxTest} I Pauli |0>", 1E-10);
            }
        }

        OptimizedBEXY(pauliBasisQubit, LittleEndian(indexRegister), targetRegister);
        Adjoint ApplyXorInPlace(targetIndex, LittleEndian(indexRegister));

        // Choose X or Y operator.
        if pauliBasis == PauliX {
            // no op
        } elif pauliBasis == PauliY {
            X(pauliBasisQubit);
        }
    }


    // Test OptimizedBEXY operator.
    operation OptimizedBEOperatorPlusTestHelper (pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int) : Unit {

        let indexRegisterSize = Ceiling(Lg(IntAsDouble(targetRegisterSize)));
        use pauliBasisQubit = Qubit();
        use indexRegister = Qubit[indexRegisterSize];
        use targetRegister = Qubit[targetRegisterSize];
        // Choose X or Y operator.
        if(pauliBasis == PauliX) {
            // no op
        } elif pauliBasis == PauliY {
            X(pauliBasisQubit);
        }

        // Create indexRegister state.
        ApplyXorInPlace(targetIndex, LittleEndian(indexRegister));

        // Initialize targetRegister states in |+>
        ApplyToEachCA(H, targetRegister);
        OptimizedBEXY(pauliBasisQubit, LittleEndian(indexRegister), targetRegister);
        for idxTest in 0..targetRegisterSize - 1 {
            let testQubit = targetRegister[idxTest];
            if(targetIndex > idxTest){
                // Test Z Pauli
                // |0> -> |0>
                // |+> -> |->
                Message($"Test Z Pauli on qubit {idxTest}");
                Diag.AssertMeasurementProbability([PauliX], [testQubit], One, 1.0, $"Error: Test {idxTest} Z Pauli |+>", 1e-10);
            }
            elif(targetIndex == idxTest){
                // Test X Pauli
                // |0> -> |1>
                // |+> -> |+>
                if(pauliBasis == PauliX){
                    Message($"Test X Pauli on qubit {idxTest}");
                    Diag.AssertMeasurementProbability([PauliX], [testQubit], Zero, 1.0, $"Error: Test {idxTest} X Pauli |+>", 1e-10);
                }

                // Test Y Pauli
                // |0> -> i|1>
                // |+> -> -i|->
                if(pauliBasis == PauliY){
                    Message($"Test Y Pauli on qubit {idxTest}");
                    Diag.AssertMeasurementProbability([PauliX], [testQubit], One, 1.0, $"Error: Test {idxTest} Y Pauli |+>", 1e-10);
                }

            }
            else{
                // Test Identitfy Pauli
                // |0> -> |0>
                // |+> -> |+>
                Message($"Test I Pauli on qubit {idxTest}");
                Diag.AssertMeasurementProbability([PauliX], [testQubit], Zero, 1.0, $"Error: Test {idxTest} I Pauli |+>", 1e-10);
            }
        }
        OptimizedBEXY(pauliBasisQubit, LittleEndian(indexRegister), targetRegister);
        ApplyToEachCA(H, targetRegister);

        (Adjoint ApplyXorInPlace)(targetIndex, LittleEndian(indexRegister));

        // Choose X or Y operator.
        if pauliBasis == PauliX {
            // no op
        } elif pauliBasis == PauliY {
            X(pauliBasisQubit);
        }
    }


    operation OptimizedBEOperatorZeroTest () : Unit {
        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;

        for idxPauli in 0 .. 1 {
            let pauliBasis = paulis[idxPauli];

            for targetIndex in 0 .. targetRegisterSize - 1 {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                OptimizedBEOperatorZeroTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }


    operation OptimizedBEOperatorPlusTest () : Unit {

        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;

        for idxPauli in 0 .. 1 {
            let pauliBasis = paulis[idxPauli];

            for targetIndex in 0 .. targetRegisterSize - 1 {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                OptimizedBEOperatorPlusTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }


    // Test phase of controlled OptimizedBEXY operator.
    operation ControlledOptimizedBEOperatorTestHelper(
        pauliBasis : Pauli, targetRegisterSize : Int, targetIndex : Int
    )
    : Unit {
        let indexRegisterSize = Ceiling(Lg(IntAsDouble(targetRegisterSize)));

        use pauliBasisQubit = Qubit();
        use indexRegister = Qubit[indexRegisterSize];
        use targetRegister = Qubit[targetRegisterSize];
        use controlQubit = Qubit();
        let testQubit = targetRegister[targetIndex];

        // Create indexRegister state.
        within {
            ApplyXorInPlace(targetIndex, LittleEndian(indexRegister));

            // Initialize control in |+> state.
            H(controlQubit);

            // Choose X or Y operator.
            if pauliBasis == PauliX {
                // Initialize testQubit state in X +1 eigenstate
                H(testQubit);
                Controlled OptimizedBEXY([controlQubit], (pauliBasisQubit, LittleEndian(indexRegister), targetRegister));
            } elif pauliBasis == PauliY {
                // Initialize testQubit state Y +1 eigenstate
                X(pauliBasisQubit);
                H(testQubit);
                S(testQubit);
                Controlled OptimizedBEXY([controlQubit], (pauliBasisQubit, LittleEndian(indexRegister), targetRegister));
            }
        } apply {
            Diag.AssertPhase(0.0, controlQubit, 1E-10);
        }
    }


    operation ControlledOptimizedBEOperatorPlusTest () : Unit {
        let paulis = [PauliX, PauliY];
        let targetRegisterSize = 7;

        for idxPauli in 0 .. 1 {
            let pauliBasis = paulis[idxPauli];

            for targetIndex in 0 .. targetRegisterSize - 1 {
                Message($"pauliBasis {pauliBasis}; targetIndex {targetIndex}");
                ControlledOptimizedBEOperatorTestHelper(pauliBasis, targetRegisterSize, targetIndex);
            }
        }
    }

    // Test SelectZ operator
    operation SelectZTest () : Unit {
        let targetRegisterSize = 7;
        let indexRegisterSize = Ceiling(Lg(IntAsDouble(targetRegisterSize)));

        use targetRegister = Qubit[targetRegisterSize];
        use indexRegister = Qubit[indexRegisterSize];
        for idxTest in 0 .. targetRegisterSize - 1 {
            H(targetRegister[idxTest]);
            ApplyXorInPlace(idxTest, LittleEndian(indexRegister));
            SelectZ(LittleEndian(indexRegister), targetRegister);
            Diag.AssertMeasurementProbability([PauliX], [targetRegister[idxTest]], One, 1.0, $"Error: Test {idxTest} X Pauli |+>", 1E-10);
            Z(targetRegister[idxTest]);
            Adjoint ApplyXorInPlace(idxTest, LittleEndian(indexRegister));
            H(targetRegister[idxTest]);
        }
    }

}


