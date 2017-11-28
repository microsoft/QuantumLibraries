// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// note that this library is still work in progress. More arithmetic to be added. 
// FIXME: complete this library 

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    
    /// # Summary 
    /// `PrepareInteger` create the quantum state corresponding to the given 
    /// integer `value` in a quantum register `qs` of `n` qubits. The encoding
    /// of the integer is the little endian encoding starting with the least 
    /// significant bit of `value` in the first qubit `qs[0]`. 
    ///
    /// # Input
    /// ## value
    /// An integer which is assumed to be non-negative. 
    /// ## qs 
    /// A quantum register which is used to store `value` in little endian encoding. 
    ///
    /// # Output
    /// An array of type `Qubit[]` that contains the little endian encoding of `value`. 
     operation PrepareInteger(value : Int, qs : Qubit[]) : () {
        body {
            let bitrepresentation = BoolArrFromPositiveInt(value, Length(qs));
            for (idx in 0..Length(qs)-1) { 
                if bitrepresentation[idx] {
                    X(qs[idx]);
                }
            }            
        }
    }

    /// # Summary 
    /// `MeasureInteger` reads out the content of a quantum register and converts 
    /// it to an integer of type `Int`. The measurement is performed with respect 
    /// to the standard computational basis, i.e., the eigenbasis of `PauliZ`.
    ///
    /// # Input
    /// ## qs 
    /// A quantum register which is assumed to be in little endian encoding. 
    ///
    /// # Output
    /// An integer of type `Int` that contains the meaured value of `qs`. 
    operation MeasureInteger(qs : Qubit[]) : Int {
        body {
            mutable results = new Result[Length(qs)]; 
            for (idx in 0..Length(qs)-1) { 
                set results[idx] = MResetZ(qs[idx]);
            }
            return PositiveIntFromResultArr(results);
        }
    }

    /// # Summary
    /// Diagonal shift circuit. This implements the diagonal unitary diag(exp(2πisk/2ⁿ):k=0,...,2ⁿ-1).
    /// # Input
    /// ## qs
    /// quantum register.
    /// ## s
    /// phase shift.
    /// FIXME: This needs to be implemented
    operation DiagonalShift (s : Int, qs : Qubit[]) : () { 
        body {
            // TODO: implement the diagonal shift circuit used in the Beauregard adder
        }
    }

    /// # Summary
    /// Given a quantum register representing an integer in the phase-domain little-endian
    /// convention by the state $\operatorname{QFT} \sum_i a_i \ket{i}$ for some complex
    /// amplitudes $\{a_i\}$, and an integer $k$, increments the state of the register to
    /// $\operatorname{QFT} \sum_i a_i \ket{i + 2^k \mod 2^n}$, where
    /// $n = \texttt{Length}(\texttt{register})$.
    ///
    /// # Input
    /// ## powerOfTwo
    /// An integer $k$ such that the phase-domain state is incremented by $2^k$.
    /// ## register
    /// A register whose state is to be incremented.
    /// FIXME: This needs to be implemented
    operation PhaseIncrementPowerOf2(powerOfTwo : Int, register : Qubit[]) : () {
        body {
            let reducedRegister = register[powerOfTwo..Length(register) - 1];
            // TODO: write this more elegantly. The +1 in idxQubit kind of kills us.
            for (idxQubit in 0..Length(reducedRegister) - 1) {
                R1Frac(1, idxQubit + 1, reducedRegister[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        adjoint controlled auto
    }
    
    /// # Summary
    /// Integer increment by a constant, based on phase rotations.
    /// This implements the unitary $\ket{x} \mapsto \ket{x + \texttt{inc}}$,
    /// where the addition on computational basis labels is performed
    /// modulo $2^N$, for $N = \texttt{Length}(\texttt{qs})$.
    ///
    /// # Input
    /// ## qs
    /// quantum register.
    /// ## inc
    /// increment.
    /// FIXME: This needs to be implemented
    operation IntegerIncrement(increment : Int, qs : Qubit[]) : () { 
        body {
            let bitrepresentation = BoolArrFromPositiveInt(increment, Length(qs));
            for (idx in 0..Length(qs)-1) { 
                if bitrepresentation[idx] {
                    X(qs[idx]);
                }
            }// TODO: implement the Beauregard adder
            // basic idea: QFT' DiagonalShift QFT = increment
        }

        adjoint auto
        controlled auto
        adjoint controlled auto
    }

    /// # Summary
    /// Integer increment by a constant, based on phase rotations.
    /// This implements the unitary $\ket{x} \mapsto \ket{x + \texttt{inc} \mod N}$.
    ///
    /// # Input
    /// ## qs
    /// quantum register.
    /// ## N
    /// modulus.
    /// ## inc
    /// increment.
    /// FIXME: This needs to be implemented
    operation ModularIncrement(N : Int, inc : Int,  qs : Qubit[]) : () {
        body {
            // TODO: implement the Beauregard modular adder
            // basic idea: https://arxiv.org/abs/quant-ph/0205095  (Section 2.2)
        }
    }
}

