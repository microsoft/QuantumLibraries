// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Logical;


    // # Summary
    /// Applies single-qubit gate defined by 2x2 unitary matrix.
    ///
    /// # Input
    /// ## u
    /// 2x2 unitary matrix describing the operation. 
    /// Assumed to be unitary. If it's not unitary, behaviour is undefined.
    /// ## qubit
    /// Qubit to which apply the operation.
    operation ApplySingleQubitUnitary(u: Complex[][], qubit : Qubit) : Unit is Adj + Ctl {  
        // ZYZ decomposition.
        let theta = ArcCos(AbsComplex(u[0][0]));
        let lmbda = ArgComplex(u[0][0]);
        let mu = ArgComplex(u[0][1]);
        if (AbsD(mu-lmbda) > 1e-9) { Rz(mu - lmbda, qubit); }
        if (AbsD(theta) > 1e-9) { Ry(-2.0 * theta, qubit); }
        if (AbsD(lmbda + mu) > 1e-9) { Rz(- lmbda - mu, qubit); }

        // If this is not special unitary, apply R1 to correct global phase.
        let det = MinusC(TimesC(u[0][0], u[1][1]), TimesC(u[0][1], u[1][0]));
        let phi = ArgComplex(det);
        if(AbsD(phi) > 1e-9) {R1(phi, qubit);}
    }

    // # Summary
    /// Applies gate defined by 2^n x 2^n unitary matrix.
    ///
    /// # Input
    /// ## unitary
    /// 2^n x 2^n unitary matrix describing the operation. 
    /// If matrix is not unitary or not of suitable size, throws an exception.
    /// ## qubit
    /// Qubits to which apply the operation - register of length n.
    operation ApplyUnitary(unitary: Complex[][], qubits : LittleEndian) : Unit is Adj + Ctl {
        ApplySingleQubitUnitary(matrix, qubits![0]);
    }    
}