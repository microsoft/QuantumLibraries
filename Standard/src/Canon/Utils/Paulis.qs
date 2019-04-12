// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Given a multi-qubit Pauli operator, applies the corresponding operation to
    /// a register.
    ///
    /// # Input
    /// ## pauli
    /// A multi-qubit Pauli operator represented as an array of single-qubit Pauli operators.
    /// ## target
    /// Register to apply the given Pauli operation on.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// ApplyPauli([PauliY, PauliZ, PauliX], target);
    /// ```
    /// and
    /// ```qsharp
    /// Y(target[0]);
    /// Z(target[1]);
    /// X(target[2]);
    /// ```
    operation ApplyPauli (pauli : Pauli[], target : Qubit[]) : Unit
    {
        body (...)
        {
            for (idxPauli in IndexRange(pauli))
            {
                let P = pauli[idxPauli];
                let targ = target[idxPauli];
                
                if (P == PauliX)
                {
                    X(targ);
                }
                elif (P == PauliY)
                {
                    Y(targ);
                }
                elif (P == PauliZ)
                {
                    Z(targ);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    
    
    /// # Summary
    /// Applies a Pauli operator on each qubit in an array if the corresponding
    /// bit of a Boolean array matches a given input.
    ///
    /// # Input
    /// ## pauli
    /// Pauli operator to apply to `qubits[idx]` where `bitsApply == bits[idx]`
    /// ## bitApply
    /// apply Pauli if bit is this value
    /// ## bits
    /// Boolean register specifying which corresponding qubit in `qubits` should be operated on
    /// ## qubits
    /// Quantum register on which to selectively apply the specified Pauli operator
    ///
    /// # Remarks
    /// The Boolean array and the quantum register must be of equal length.
    operation ApplyPauliFromBitString (pauli : Pauli, bitApply : Bool, bits : Bool[], qubits : Qubit[]) : Unit
    {
        body (...)
        {
            let nBits = Length(bits);
            
            //FailOn (nbits != Length(qubits), "Number of control bits must be equal to number of control qubits")
            for (idxBit in 0 .. nBits - 1)
            {
                if (bits[idxBit] == bitApply)
                {
                    ApplyPauli([pauli], [qubits[idxBit]]);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    
    
    /// # Summary
    /// Given a single-qubit Pauli operator and the index of a qubit,
    /// returns a multi-qubit Pauli operator with the given single-qubit
    /// operator at that index and `PauliI` at every other index.
    ///
    /// # Input
    /// ## pauli
    /// A single-qubit Pauli operator to be placed at the given location.
    /// ## location
    /// An index such that `output[location] == pauli`, where `output` is
    /// the output of this function.
    /// ## n
    /// Length of the array to be returned.
    ///
    /// # Remarks
    /// ## Example
    /// To obtain the array `[PauliI, PauliI, PauliX, PauliI]`:
    /// ```qsharp
    /// EmbedPauli(PauliX, 2, 3);
    /// ```
    function EmbedPauli (pauli : Pauli, location : Int, n : Int) : Pauli[]
    {
        mutable pauliArray = new Pauli[n];
        
        for (index in 0 .. n - 1)
        {
            if (index == location)
            {
                set pauliArray[index] = pauli;
            }
            else
            {
                set pauliArray[index] = PauliI;
            }
        }
        
        return pauliArray;
    }
    
    
    // FIXME: Remove in favor of something that computes arbitrary
    //        weight Paulis.
    
    /// # Summary
    /// Returns an array of all weight-1 Pauli operators
    /// on a given number of qubits.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits on which the returned Pauli operators
    /// are defined.
    ///
    /// # Output
    /// An array of multi-qubit Pauli operators, each of which is
    /// represented as an array with length `nQubits`.
    function WeightOnePaulis (nQubits : Int) : Pauli[][]
    {
        mutable paulis = new Pauli[][3 * nQubits];
        let pauliGroup = [PauliX, PauliY, PauliZ];
        
        for (idxQubit in 0 .. nQubits - 1)
        {
            for (idxPauli in IndexRange(pauliGroup))
            {
                set paulis[idxQubit * Length(pauliGroup) + idxPauli] = EmbedPauli(pauliGroup[idxPauli], idxQubit, nQubits);
            }
        }
        
        return paulis;
    }

}


