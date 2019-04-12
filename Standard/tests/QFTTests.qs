// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    // To test QFT we hard code circuits based on Figure 5.1 on Page 219 of
    // [ *Michael A. Nielsen , Isaac L. Chuang*,
    //    Quantum Computation and Quantum Information ](http://doi.org/10.1017/CBO9780511976667)

    /// # Summary
    /// Hard-code 1 qubit QFT
    operation QFT1 (target : BigEndian) : Unit {
        body (...) {
            EqualityFactI(Length(target!), 1, $"`Length(target!)` must be 1");
            H((target!)[0]);
        }

        adjoint invert;
    }

    /// # Summary
    /// Hard-code 2 qubit QFT
    operation QFT2 (target : BigEndian) : Unit {
        body (...) {
            EqualityFactI(Length(target!), 2, $"`Length(target!)` must be 2");
            let (q1, q2) = ((target!)[0], (target!)[1]);
            H(q1);
            Controlled R1Frac([q2], (2, 2, q1));
            H(q2);
            SWAP(q1, q2);
        }

        adjoint invert;
    }

    /// # Summary
    /// Hard-code 3 qubit QFT
    operation QFT3 (target : BigEndian) : Unit {
        body (...) {
            EqualityFactI(Length(target!), 3, $"`Length(target)` must be 3");
            let (q1, q2, q3) = ((target!)[0], (target!)[1], (target!)[2]);
            H(q1);
            Controlled R1Frac([q2], (2, 2, q1));
            Controlled R1Frac([q3], (2, 3, q1));
            H(q2);
            Controlled R1Frac([q3], (2, 2, q2));
            H(q3);
            SWAP(q1, q3);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Hard-code 4 qubit QFT
    operation QFT4 (target : BigEndian) : Unit {
        
        body (...) {
            EqualityFactI(Length(target!), 4, $"`Length(target!)` must be 4");
            let (q1, q2, q3, q4) = ((target!)[0], (target!)[1], (target!)[2], (target!)[3]);
            H(q1);
            Controlled R1Frac([q2], (2, 2, q1));
            Controlled R1Frac([q3], (2, 3, q1));
            Controlled R1Frac([q4], (2, 4, q1));
            H(q2);
            Controlled R1Frac([q3], (2, 2, q2));
            Controlled R1Frac([q4], (2, 3, q2));
            H(q3);
            Controlled R1Frac([q4], (2, 2, q3));
            H(q4);
            SWAP(q1, q4);
            SWAP(q2, q3);
        }
        
        adjoint invert;
    }
    
    
    operation ApplyBEToRegisterA (op : (BigEndian => Unit : Adjoint), target : Qubit[]) : Unit {
        
        body (...) {
            op(BigEndian(target));
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Compares QFT to the hard-coded implementations
    operation QFTTest () : Unit {
        let testFunctions = [QFT1, QFT2, QFT3, QFT4];

        for (i in IndexRange(testFunctions)) {
            AssertOperationsEqualReferenced(i + 1, ApplyBEToRegisterA(testFunctions[i], _), ApplyBEToRegisterA(QFT, _));
        }
    }

}

