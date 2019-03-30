// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    open Microsoft.Quantum.Primitive;
    
    
    /// # Summary
    /// Measures the given Pauli operator using an explicit scratch
    /// qubit to perform the measurement.
    ///
    /// # Input
    /// ## pauli
    /// A multi-qubit Pauli operator specified as an array of
    /// single-qubit Pauli operators.
    /// ## target
    /// Qubit register to be measured.
    ///
    /// # Output
    /// The result of measuring the given Pauli operator on
    /// the `target` register.
    operation MeasureWithScratch (pauli : Pauli[], target : Qubit[]) : Result
    {
        mutable result = Zero;
        
        using (scratchRegister = Qubit[1])
        {
            let scratch = scratchRegister[0];
            H(scratch);
            
            for (idxPauli in 0 .. Length(pauli) - 1)
            {
                let P = pauli[idxPauli];
                let src = target[idxPauli];
                
                if (P == PauliX)
                {
                    Controlled X([scratch], src);
                }
                elif (P == PauliY)
                {
                    Controlled Y([scratch], src);
                }
                elif (P == PauliZ)
                {
                    Controlled Z([scratch], src);
                }
            }
            
            H(scratch);
            set result = M(scratch);
            ResetAll(scratchRegister);
        }
        
        return result;
    }
    
    
    /// # Summary
    /// Returns one of the single-qubit Pauli operators uniformly
    /// at random.
    ///
    /// # Output
    /// A `Pauli` operator that is one of `[PauliI, PauliX, PauliY, PauliZ]`.
    ///
    /// # Remarks
    /// This function calls <xref:microsoft.quantum.primitive.random>, so
    /// its randomness depends on the implementation of `Random`.
    operation RandomSingleQubitPauli () : Pauli
    {
        let probs = [0.5, 0.5, 0.5, 0.5];
        let idxPauli = Random(probs);
        let singleQubitPaulis = [PauliI, PauliX, PauliY, PauliZ];
        return singleQubitPaulis[idxPauli];
    }
    
    
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
            for (idxPauli in 0 .. Length(pauli) - 1)
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
    /// Given an array of multi-qubit Pauli operators, measures each using a specified measurement
    /// gadget, then returns the array of results.
    ///
    /// # Input
    /// ## paulis
    /// Array of multi-qubit Pauli operators to measure.
    /// ## target
    /// Register on which to measure the given operators.
    /// ## gadget
    /// Operation which performs the measurement of a given multi-qubit operator.
    ///
    /// # Output
    /// The array of results obtained from measuring each element of `paulis`
    /// on `target`.
    operation MeasurePaulis (paulis : Pauli[][], target : Qubit[], gadget : ((Pauli[], Qubit[]) => Result)) : Result[]
    {
        mutable results = new Result[Length(paulis)];
        
        for (idxPauli in 0 .. Length(paulis) - 1)
        {
            set results[idxPauli] = gadget(paulis[idxPauli], target);
        }
        
        return results;
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
            for (idxPauli in 0 .. Length(pauliGroup) - 1)
            {
                set paulis[idxQubit * Length(pauliGroup) + idxPauli] = EmbedPauli(pauliGroup[idxPauli], idxQubit, nQubits);
            }
        }
        
        return paulis;
    }
    
    
    // NB: This operation is intended to be private to Paulis.qs.
    
    operation _BasisChangeZtoY (target : Qubit) : Unit
    {
        body (...)
        {
            H(target);
            S(target);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Measures each qubit in a given array in the standard basis.
    /// # Input
    /// ## targets
    /// An array of qubits to be measured.
    /// # Output
    /// An array of measurement results.
    operation MultiM (targets : Qubit[]) : Result[]
    {
        mutable results = new Result[Length(targets)];
        
        for (idxQubit in 0 .. Length(targets) - 1)
        {
            set results[idxQubit] = M(targets[idxQubit]);
        }
        
        return results;
    }
    
    
    /// # Summary
    /// Measures a single qubit in the `Z` basis,
    /// and resets it to the standard basis state
    /// |0〉 following the measurement.
    ///
    /// # Input
    /// ## target
    /// A single qubit to be measured.
    ///
    /// # Output
    /// The result of measuring `target` in the Pauli $Z$ basis.
    operation MResetZ (target : Qubit) : Result
    {
        let result = M(target);
        
        if (result == One)
        {
            // Recall that the +1 eigenspace of a measurement operator corresponds to
            // the Result case Zero. Thus, if we see a One case, we must reset the state
            // have +1 eigenvalue.
            X(target);
        }
        
        return result;
    }
    
    
    /// # Summary
    /// Measures a single qubit in the X basis,
    /// and resets it to the standard basis state
    /// |0〉 following the measurement.
    ///
    /// # Input
    /// ## target
    /// A single qubit to be measured.
    ///
    /// # Output
    /// The result of measuring `target` in the Pauli $X$ basis.
    operation MResetX (target : Qubit) : Result
    {
        let result = Measure([PauliX], [target]);
        
        // We must return the qubit to the Z basis as well.
        H(target);
        
        if (result == One)
        {
            // Recall that the +1 eigenspace of a measurement operator corresponds to
            // the Result case Zero. Thus, if we see a One case, we must reset the state
            // have +1 eigenvalue.
            X(target);
        }
        
        return result;
    }
    
    
    /// # Summary
    /// Measures a single qubit in the Y basis,
    /// and resets it to the standard basis state
    /// |0〉 following the measurement.
    ///
    /// # Input
    /// ## target
    /// A single qubit to be measured.
    ///
    /// # Output
    /// The result of measuring `target` in the Pauli $Y$ basis.
    operation MResetY (target : Qubit) : Result
    {
        let result = Measure([PauliY], [target]);
        
        // We must return the qubit to the Z basis as well.
        Adjoint _BasisChangeZtoY(target);
        
        if (result == One)
        {
            // Recall that the +1 eigenspace of a measurement operator corresponds to
            // the Result case Zero. Thus, if we see a One case, we must reset the state
            // have +1 eigenvalue.
            X(target);
        }
        
        return result;
    }
    
    
    /// # Summary
	/// Applies the Y-basis analog to the Hadamard transformation
	/// that interchanges the Z and Y axes.
	/// 
    /// The Y Hadamard transformation $H_Y = S H$ on a single qubit is:
    ///
    /// \begin{align}
    ///     H_Y \mathrel{:=}
    ///     \frac{1}{\sqrt{2}}
    ///     \begin{bmatrix}
    ///         1 & 1 \\\\
    ///         i & -i
    ///     \end{bmatrix}.
    /// \end{align}
    ///
    /// # Input
    /// ## qubit
    /// Qubit to which the gate should be applied.
    ///
    /// # See Also
    operation HY (target : Qubit) : Unit
    {
        body (...)
        {
            H(target);
            S(target);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
}


