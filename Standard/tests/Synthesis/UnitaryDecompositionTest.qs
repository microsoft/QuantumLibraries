// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Synthesis;
    open Microsoft.Quantum.Random;
    
    operation ApplyUnitaryToRegister(matrix: Complex[][], qubits: Qubit[]) : Unit {
        ApplyUnitary(matrix, LittleEndian(qubits));
    }

    operation CheckOperation(matrix: Complex[][], expected: (Qubit[] => Unit is Adj)) : Unit {
        let nQubits = 1; // TODO: number of qubits is not always 1.
        AssertOperationsEqualInPlace(nQubits, ApplyUnitaryToRegister(matrix, _), expected);
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_Identity () : Unit {
        let matrix = [[Complex(1.0, 0.0), Complex(0.0, 0.0)],
                      [Complex(0.0, 0.0), Complex(1.0, 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(I, _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_PauliX () : Unit {
        let matrix = [[Complex(0.0, 0.0), Complex(1.0, 0.0)],
                      [Complex(1.0, 0.0), Complex(0.0, 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(X, _));
    }


    @Test("QuantumSimulator")
    operation ApplyUnitary_PauliY () : Unit {
        let matrix = [[Complex(0.0, 0.0), Complex(0.0, -1.0)],
                      [Complex(0.0, 1.0), Complex(0.0, 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(Y, _));
    }


    @Test("QuantumSimulator")
    operation ApplyUnitary_PauliZ () : Unit {
        let matrix = [[Complex(1.0, 0.0), Complex(0.0, 0.0)],
                      [Complex(0.0, 0.0), Complex(-1.0, 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(Z, _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_Hadamard () : Unit {
        let matrix = [[Complex(Sqrt(0.5), 0.0), Complex(Sqrt(0.5), 0.0)],
                      [Complex(Sqrt(0.5), 0.0), Complex(-Sqrt(0.5), 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(H, _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_HadamardY () : Unit {
        let matrix = [[Complex(Sqrt(0.5), 0.0), Complex(Sqrt(0.5), 0.0)],
                      [Complex(0.0, Sqrt(0.5)), Complex(0.0, -Sqrt(0.5))]];
        CheckOperation(matrix, ApplyToHeadA(HY, _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_Rx () : Unit {
        let matrix = [[Complex(Cos(1.0), 0.0), Complex(0.0, -Sin(1.0))],
                      [Complex(0.0, -Sin(1.0)), Complex(Cos(1.0), 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(Rx(2.0, _), _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_Ry () : Unit {
        let matrix = [[Complex(Cos(1.0), 0.0), Complex(-Sin(1.0), 0.0)],
                      [Complex(Sin(1.0), 0.0), Complex(Cos(1.0), 0.0)]];
        CheckOperation(matrix, ApplyToHeadA(Ry(2.0, _), _));
    }
}
