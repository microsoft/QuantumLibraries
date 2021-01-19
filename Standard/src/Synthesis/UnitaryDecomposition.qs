// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Applies single-qubit gate defined by 2x2 unitary matrix.
    internal operation ApplySingleQubitUnitary(u : Complex[][], qubit : Qubit) : Unit is Adj + Ctl {  
        // ZYZ decomposition.
        let theta = ArcCos(AbsComplex(u[0][0]));
        let lmbda = ArgComplex(u[0][0]);
        let mu = ArgComplex(u[0][1]);
        if (AbsD(mu - lmbda) > 1e-10) { Rz(mu - lmbda, qubit); }
        if (AbsD(theta) > 1e-10) { Ry(-2.0 * theta, qubit); }
        if (AbsD(lmbda + mu) > 1e-10) { Rz(-lmbda - mu, qubit); }

        // If this is not special unitary, apply R1 to correct the global phase.
        let det = MinusC(TimesC(u[0][0], u[1][1]), TimesC(u[0][1], u[1][0]));
        let phi = ArgComplex(det);
        if (AbsD(phi) > 1e-10) { R1(phi, qubit); }
    }

    function _TwoLevelDecomposition(unitary: Complex[][]) : (Complex[][], Int, Int)[] {
        body intrinsic;
    }

    /// # Summary
    /// For every 2-level unitary calculates "flip mask", which denotes qubits which should 
    /// be inverted before and after applying corresponding 1-qubit gate.
    /// For convenience prepends result with 0.
    internal function FlipMasks(decomposition: (Complex[][], Int, Int)[], 
                                allQubitsMask: Int) : Int[] {
        let n = Length(decomposition);
        mutable flipMasks = ConstantArray(n + 1, 0);
        for ((i, (_, i1, i2)) in Enumerated(decomposition)) {
            set flipMasks w/= (i + 1) <- (allQubitsMask - i2); 
        }
        return flipMasks;
    }

    /// # Summary
    /// Applies gate defined by 2^n x 2^n unitary matrix.
    ///
    /// Fails if matrix is not unitary, or has wrong size. 
    ///
    /// # Input
    /// ## unitary
    /// 2^n x 2^n unitary matrix describing the operation. 
    /// If matrix is not unitary or not of suitable size, throws an exception.
    /// ## qubits
    /// Qubits to which apply the operation - register of length n.
    operation ApplyUnitary(unitary: Complex[][], qubits : LittleEndian) : Unit is Adj + Ctl {
        SquareMatrixFact(unitary);
        EqualityFactI(Length(unitary), 1 <<< Length(qubits!),
            "Matrix size is not consistent with register length.");

        let allQubitsMask = (1 <<< Length(qubits!)) - 1;
        let decomposition = _TwoLevelDecomposition(unitary);
        let flipMasks = FlipMasks(decomposition, allQubitsMask);
        
        // i1, i2 - indices of non-trivial 2x2 submatrix of two-level unitary matrix being 
        // applied. 
        // i1 and i2 differ in exactly one bit; i1 < i2.
        // matrix - 2x2 non-trivial unitary submatrix of said two-level unitary.
        for ((i, (matrix, i1, i2)) in Enumerated(decomposition)) {
            ApplyXorInPlace(flipMasks[i + 1] ^^^ flipMasks[i], qubits);

            let targetMask = i1 ^^^ i2;
            let controlMask = allQubitsMask - targetMask;
            let (controls, targets) = MaskToQubitsPair(qubits!, MCMTMask(controlMask, targetMask));
            Controlled ApplySingleQubitUnitary(controls, (matrix, targets[0])); 
        }   
        ApplyXorInPlace(Tail(flipMasks), qubits);
    }    

    /// # Summary
    /// Checks that given array represents a square matrix.
    internal function SquareMatrixFact(matrix : Complex[][]) : Unit {
        let n = Length(matrix);
        for (row in matrix) {
            EqualityFactI(Length(row), n, "Matrix is not square.");
        }
    } 
}
