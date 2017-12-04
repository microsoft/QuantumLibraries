// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Shor {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    //////////////////////////////////////////////////////////////////////////
    // Introduction //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    /// This sample contains Q# code implementing Shor's quantum algorithm for
    /// factoring integers. The underlying modular arithmetic is implemented 
    /// in phase encoding, based on paper by Stephane Beauregard who gave a
    /// quantum circuit for factoring n-bit numbers that needs 2n+3 qubits and 
    /// O(n^3 log(n)) many elementary quantum gates. The Q# implementation 
    /// shown here is closely patterned after a implementation in the 
    /// LIQUi|> quantum computing programming language. 

    //////////////////////////////////////////////////////////////////////////
    // Util functions to create and measure quantum integer states ///////////
    //////////////////////////////////////////////////////////////////////////
    
    /// # Summary 
    /// `R` is a single qubit unitary operation that applies a diagonal operation to
    /// the state of the input qubit `q` around the `PauliZ` axis. Mathematically, the 
    /// operator $R$ corresponds to the unitary $R = {\rm diag}(1,\exp(2\pi i/2^k)$. 
    ///
    /// # Input
    /// ## k 
    /// An integer value which is used to specify the applied diagonal operator. 
    /// ## q 
    /// A qubit which is rotated by $R = {\rm diag}(1,\exp(2\pi i/2^k)$.
    ///
    /// # Output
    /// The state after the rotation $R$ has been performed. 
    operation R(k : Int, q : Qubit) : () {
        body { 
            R1Frac(2, k, q);
        }
        adjoint auto
        controlled auto
        adjoint controlled auto
    }

    /// # Summary 
    /// `CR` implements a single controlled qubit unitary operation that applies a diagonal operation 
    /// to the state of two input qubits `c` and `t`. Mathematically, the operator $CR$ corresponds 
    /// to the unitary $CR = {\rm diag}(1, 1, 1, \exp(2\pi i/2^k)$. 
    ///
    /// # Input
    /// ## k 
    /// An integer value which is used to specify the element $(3,3)$ of the applied diagonal operator. 
    /// ## c
    /// A qubit which serves as the control qubit of `CR`. 
    /// ## t
    /// A qubit which serves as the target qubit of `CR`. 
    ///
    /// # Output
    /// The state after the controlled rotation $CR$ has been performed. 
    operation CR(k : Int, c : Qubit, t : Qubit) : () {
        body { 
            (Controlled R(k, _))([c], t); 
        }
        adjoint auto
    }

    /// # Summary 
    /// `CCR` implements a doubly controlled qubit unitary operation that applies a diagonal operation 
    /// to the state of three input qubits `c1`, `c2` and `t`. Mathematically, the operator $CCR$ corresponds 
    /// to the $8 \times 8$ unitary $CCR = {\rm diag}(1, 1, \ldots, 1, \exp(2\pi i/2^k)$. 
    ///
    /// # Input
    /// ## k 
    /// An integer value which is used to specify the element $(7,7)$ of the applied diagonal operator. 
    /// ## c1
    /// A qubit which serves as the first control qubit of `CR`. 
    /// ## c2
    /// A qubit which serves as the second control qubit of `CR`. 
    /// ## t
    /// A qubit which serves as the target qubit of `CR`. 
    ///
    /// # Output
    /// The state after the controlled rotation $CCR$ has been performed. 
    operation CCR(k : Int, c1 : Qubit, c2 : Qubit, t : Qubit) : () {
        body {
            (Controlled R(k, _))([c1; c2], t);
        }
        adjoint auto
    }

    /// # Summary 
    /// An operation to add a classical constant to a quantum register. The quantum register
    /// is assumed to be in phase encoding, i.e., a quantum Fourier transform (QFT) was 
    /// applied to an integer that itself was encoded in little endian encoding. 
    /// Mathematically, the operation corresponds to $\ket{x} \mapto \ket{x+a}$
    /// Note that this operation can also perform the subtraction of a value from a 
    /// quantum register when run in reverse, i.e., when applying `(Adjoint opGen)`. 
    ///
    /// # Input
    /// ## a
    /// A constant integer value which is used to specify the increment. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    /// ## ctrlCnt
    /// An integer that specifies whether the addition is performed regularly (`ctrlCnt=0`)
    /// or controlled on 1 control qubit (`ctrlCnt=0`) or controlled on 2 control qubits 
    /// (`ctrlCnt=2`). 
    ///
    /// # Output
    /// The result of adding the constant `a` to the integer `qs` in phase encoding. 
    operation opGen(a : Int, qs : Qubit[], ctrlCnt : Int) : () { 
        body { 
            let n   = Length(qs)-ctrlCnt;
            let cs  = qs[0..ctrlCnt-1];
            let bs  = qs[ctrlCnt..Length(qs)-1];
            let aBitRep = BoolArrFromPositiveInt(a,n);
            
            mutable k = 0; 
            for (bIdx in 0..n-1) { 
                for (aIdx in 0..bIdx) {
                    if aBitRep[aIdx] {
                        set k = 1+bIdx-aIdx;
                        if (ctrlCnt == 0) {
                            R(k, bs[bIdx]); 
                        }
                        elif (ctrlCnt == 1) {
                            CR(k, cs[0], bs[bIdx]);
                        }
                        else { 
                            CCR(k, cs[0], cs[1], bs[bIdx]);
                        }
                    }
                }
            }
        }
        adjoint auto
    }

    /// # Summary 
    /// An operation to add a classical constant to a quantum register in phase encoding. 
    ///
    /// # Input
    /// ## a
    /// A constant integer value which is used to specify the increment. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    ///
    /// # Output
    /// The result of adding the constant `a` to the integer `qs` in phase encoding. 
    operation AddA(a : Int, qs : Qubit[]) : () {
        body { 
            opGen(a, qs, 0); 
        }
        adjoint auto
    }
    
    /// # Summary 
    /// Controlled addition of a classical constant to a quantum register in phase encoding. 
    ///
    /// # Input
    /// ## a
    /// A constant integer value which is used to specify the increment. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. The first qubit 
    /// of `qs` is assumed to be the control qubit and the remaining qubits of `qs` are
    /// assumed to represent the integer in little endian phase encoding.
    ///
    /// # Output
    /// The result of controlled addition of constant `a` to integer `qs` in phase encoding. 
    operation CAddA(a : Int, qs : Qubit[]) : () { 
        body { 
            opGen(a, qs, 1); 
        }
        adjoint auto
    }
    
    /// # Summary 
    /// Doubly-controlled addition of a classical constant to a quantum register in phase encoding. 
    ///
    /// # Input
    /// ## a
    /// A constant integer value which is used to specify the increment. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. The first two qubits
    /// of `qs` is assumed to be the control qubits and the remaining qubits of `qs` are 
    /// assumed to represent the integer in little endian phase encoding.
    ///
    /// # Output
    /// The result of doubly-controlled addition of constant `a` to integer `qs` in phase encoding. 
    operation CCAdd(a : Int, qs : Qubit[]) : () { 
        body { 
            opGen(a, qs, 2); 
        }
        adjoint auto
    }

    /// # Summary 
    /// Quantum circuit to perform a quantum Fourier transform (QFT) to a quantum register. 
    /// It should be noted that this operation slightly differs from applying from the canon
    /// opeation `QFT`. When applied to a quantum register `qs` the difference is that `QFT` 
    /// performs the final bit-reversal permutation whereas `QFTcirc` does not. 
    ///
    /// # Input
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    ///
    /// # Output
    /// The result of applying the quantum Fourier transform (sans bit-reversal) to `qs`. 
    operation QFTcirc(qs : Qubit[]) : () { 
        body { 
            let n = Length(qs);
            for (aIdx in n-1..-1..0) { 
                H(qs[aIdx]);
                for (k in 2..aIdx+1) { 
                    let c = qs[aIdx-(k-1)];
                    CR(k, c, qs[aIdx]);
                }
            }
        }
        adjoint auto
    }
    
    /// # Summary 
    /// An operation to add a classical constant to a quantum register where the 
    /// underlying numbers are reduced modulo a given classical constant. The quantum 
    /// register is assumed to be in phase encoding, i.e., a quantum Fourier transform 
    /// (QFT) was applied to an integer that itself was encoded in little endian encoding. 
    /// Mathematically, the operation corresponds to $\ket{x} \mapto \ket{(x+a) \% m}$, 
    /// where $m$ is the modulus, specified here as the input parameter `modulus`.
    /// Note that this operation can also perform the subtraction of a value from a 
    /// quantum register by running the operation in revese, i.e., by applying the
    /// operation `(Adjoint AddModN)`.
    ///
    /// # Input
    /// ## modulus
    /// A constant integer value which is used to specify the modulus with respect to 
    /// which the increment by constant `a` is performed.
    /// ## a
    /// A constant integer value which is used to specify the increment. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    ///
    /// # Output
    /// The result of adding the constant `a` to the integer `qs` in phase encoding,
    /// where the operation is performed modulo the integer `modulus`.
    operation AddModN(modulus : Int, a : Int, qs : Qubit[]) : () { 
        body { 
            let n = Length(qs)-4;
            let c1 = qs[0];
            let c2 = qs[1];
            let bs = qs[2..n+2];
            let anc = qs[Length(qs)-1];
            let bMx = bs[Length(bs)-1];
            let cbs = [c1;c2] + bs; 
        
            CCAdd(a, cbs); 
            (Adjoint AddA)(modulus, bs); 
            (Adjoint QFTcirc)(bs);
            CNOT(bMx, anc); 
            QFTcirc(bs); 
            CAddA(modulus, [anc]+bs); 
            (Adjoint CCAdd)(a, cbs);
            (Adjoint QFTcirc)(bs); 
            X(bMx);
            CNOT(bMx, anc);
            X(bMx); 
            QFTcirc(bs); 
            CCAdd(a, cbs);
        }
        adjoint auto
    }

    /// # Summary 
    /// An operation to multiply a classical constant with the number in a quantum register,
    /// where the underlying numbers are reduced modulo a given classical constant. The quantum 
    /// register is assumed to be in phase encoding, i.e., a quantum Fourier transform 
    /// (QFT) was applied to an integer that itself was encoded in little endian encoding. 
    /// Mathematically, the operation corresponds to $\ket{x} \mapto \ket{(a x) \% m}$, 
    /// where $m$ is the modulus, specified here as the input parameter `modulus`.
    ///
    /// # Input
    /// ## modulus
    /// A constant integer value which is used to specify the modulus with respect to 
    /// which the multiplication by constant `a` is performed.
    /// ## a
    /// A constant integer value which is used to specify the constant multiplication. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    ///
    /// # Output
    /// The result of adding the constant `a` to the integer `qs` in phase encoding,
    /// where the multiplication operation is performed modulo the integer `modulus`.
    operation MulModN(modulus : Int, a : Int, qs : Qubit[]) : () { 
        body { 
            let n = (Length(qs) - 3) / 2;
            let c = qs[0];
            let xs = qs[1..n];
            let bs = qs[(n+1)..(2*n+1)];
            let anc = qs[Length(qs)-1];
            
            mutable a2 = 0; 
            QFTcirc(bs);
            for (j in 0..n-1) {
                set a2 = (2^j * a) % modulus;
                let args = [c] + [xs[j]] + bs + [anc];
                AddModN(modulus, a2, args);
            }
            (Adjoint QFTcirc)(bs);
        }
        adjoint auto
    }

    /// # Summary 
    /// Controlled swap operation of the contents of two qubits. 
    ///
    /// # Input
    /// ## control
    /// A control qubit. The swap operation of `a` and `b` is performed if and only if 
    /// the content of `control` is in the $\ket{0}$ state. 
    /// ## a
    /// An object of type `Qubit` which is swapped with another object `b` of the same type. 
    /// ## b
    /// An object of type `Qubit` which is swapped with another object `a` of the same type. 
    ///
    /// # Output
    /// The result of applying the controlled swap operation of `a` and `b`.
    operation CSWAP(control : Qubit, a : Qubit, b : Qubit) : () { 
        body { 
            CNOT(b, a); 
            CCNOT(control, a, b); 
            CNOT(b, a); 
        }
    }    

    /// # Summary 
    /// A function used to perform one step in the modular exponentiation of a quantum 
    /// register. The encoding of the integers in the quantum register is little 
    /// endian based phase encoding, the computation is performed modulo a classical 
    /// constant value `modulus` and the constant that is part of the exponentiation
    /// is given by `a`. 
    ///
    /// # Input
    /// ## modulus
    /// A constant integer value which is used to specify the modulus with respect to 
    /// which the increment by constant `a` is performed.
    /// ## a
    /// A constant integer value which is used to specify the constant multiplication. 
    /// ## aRec
    /// A constant integer value which is used to specify the constant multiplication
    /// by $a^{-1}$, where this computation is performed modulo the given `modulus`. 
    /// The value of `aRec` can be precomputed in a classical precomputation that 
    /// happens outside the Q# framework. 
    /// ## qs
    /// An array of qubits that represents an integer in phase encoding. 
    /// # Output
    /// The result of applying the controlled multiplication of `a` to the register. 
    operation Ua(modulus : Int, a : Int, aRec : Int, qs : Qubit[]) : () { 
        body {
            let n = (Length(qs) - 3) / 2;
            let c = qs[0];
            let xs = qs[1..n];
            let bs = qs[(n+1)..(2*n+1)];
            MulModN(modulus, a, qs);
            for (i in 0..n-1) { 
                CSWAP(c, xs[i], bs[i]);
            }
            (Adjoint MulModN)(modulus, aRec, qs); 
        }
    }
}

