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
    
    internal operation ApplyUnitaryToRegister(matrix: Complex[][], qubits: Qubit[]) : Unit {
        ApplyUnitary(matrix, LittleEndian(qubits));
    }

    // Checks that `ApplyUnitary(matrix)` is equvalent to `expected` operation.
    internal operation CheckOperation(matrix: Complex[][], expected: (Qubit[] => Unit is Adj)) : Unit {
        let nQubits = Floor(Lg(IntAsDouble(Length(matrix))));
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

    @Test("QuantumSimulator")
    operation ApplyUnitary_CNOT () : Unit {
        let ZERO = Complex(0.0, 0.0);
        let ONE = Complex(1.0, 0.0);
        // Matrix for CNOT(q[0], q[1]).
        let matrix = [[ONE, ZERO, ZERO, ZERO],
                      [ZERO, ZERO, ZERO, ONE],
                      [ZERO, ZERO, ONE, ZERO],
                      [ZERO, ONE, ZERO, ZERO]];
        CheckOperation(matrix, ApplyToFirstTwoQubitsA(CNOT, _));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_SWAP () : Unit {
        let ZERO = Complex(0.0, 0.0);
        let ONE = Complex(1.0, 0.0);
        let matrix = [[ONE, ZERO, ZERO, ZERO],
                      [ZERO, ZERO, ONE, ZERO],
                      [ZERO, ONE, ZERO, ZERO],
                      [ZERO, ZERO, ZERO, ONE]];
        CheckOperation(matrix, ApplyToFirstTwoQubitsA(SWAP, _));
    }

    internal operation ApplyQFT(qubits: Qubit[]) : Unit is Adj {
        QFT(BigEndian(Reversed(qubits)));
    }

    @Test("QuantumSimulator")
    operation ApplyUnitary_QFT () : Unit {
        let matrix = [[Complex(0.5, 0.0), Complex(0.5, 0.0), Complex(0.5, 0.0), Complex(0.5, 0.0)],
                      [Complex(0.5, 0.0), Complex(0.0, 0.5), Complex(-0.5, 0.0), Complex(0.0, -0.5)],
                      [Complex(0.5, 0.0), Complex(-0.5, 0.0), Complex(0.5, 0.0), Complex(-0.5, 0.0)],
                      [Complex(0.5, 0.0), Complex(0.0, -0.5), Complex(-0.5, 0.0), Complex(0.0, 0.5)]];
        CheckOperation(matrix, ApplyQFT);
    }
    
    @Test("QuantumSimulator")
    operation ApplyUnitary_CCNOT () : Unit {
        // Matrix for CCNOT(q[0], q[1], q[2]).
        let ZERO = Complex(0.0, 0.0);
        let ONE = Complex(1.0, 0.0);
        let matrix = [[ONE, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO],
                      [ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO],
                      [ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO],
                      [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE],
                      [ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO],
                      [ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO],
                      [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO],
                      [ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO]];
        CheckOperation(matrix, ApplyToFirstThreeQubitsA(CCNOT, _));
    }
}
