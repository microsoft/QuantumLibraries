// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Chemistry;

	/// # Summary
    /// Computes Z component of Jordan-Wigner string between
    /// fermion indices in a fermionic operator with an even
    /// number of creation / annihilation operators.
    ///
    /// # Input
    /// ## nFermions
    /// The Number of fermions in the system.
    /// ## idxFermions
    /// fermionic operator indices.
    ///
    /// # Output
    /// Bitstring `Bool[]` that is `true` where a `PauliZ` should be applied.
    ///
    /// # Example
    /// let bitString = _ComputeJordanWignerBitString(6, [0,1,2,6]) ;
    /// // bitString is [false, false, false ,true, true, true, false].
    function _ComputeJordanWignerBitString(nFermions: Int,  idxFermions: Int[]) : Bool[] {
        if(Length(idxFermions) % 2 != 0){
            fail $"ComputeJordanWignerString failed. `idxFermions` must contain an even number of terms.";
        }
        
        mutable zString = new Bool[nFermions];
        for (fermionIdx in idxFermions){
            if(fermionIdx >= nFermions){
                fail $"ComputeJordanWignerString failed. fermionIdx {fermionIdx} out of range.";
        }
        for(idx in 0..fermionIdx){
            set zString[idx] = zString[idx] ? false | true;
        }
        }
        
        for (fermionIdx in idxFermions){
            set zString[fermionIdx] = false;
        }
        return zString;
    }

    // Identical to `_ComputeJordanWignerBitString`, except with the map  
    // false -> PauliI and true -> PauliZ
    function _ComputeJordanWignerPauliZString(nFermions: Int,  idxFermions: Int[]) : Pauli[] {
        let bitString = _ComputeJordanWignerBitString(nFermions, idxFermions);
        return PauliFromBitString (PauliZ, true, bitString);
    }

    // Identical to `_ComputeJordanWignerPauliZString`, except that some
    // specified elements are substituted.
    function _ComputeJordanWignerPauliString(nFermions: Int,  idxFermions: Int[], pauliReplacements : Pauli[]) : Pauli[] {
        mutable pauliString = _ComputeJordanWignerPauliZString(nFermions, idxFermions);

        for(idx in 0..Length(idxFermions)-1){
            let idxFermion = idxFermions[idx];
            let op = pauliReplacements[idx];
            set pauliString[idxFermion] = op;
        }
        
        return pauliString;
    }

    /// # Summary
    /// Applies time-evolution by a cluster operator PQ term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a cluster operator PQ term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerClusterOperatorPQTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let p = idxFermions[0];
            let q = idxFermions[1];
            let angle = 0.5 * coeff[0] * stepSize;
            let ops = [[PauliX, PauliY], [PauliY, PauliX]];
            let signs = [+1.0, -1.0];

            for (idxOp in 0 .. Length(ops) - 1) {
                let pauliString = _ComputeJordanWignerPauliString(Length(qubits), idxFermions, ops[idxOp]);
                let sign = signs[idxOp];
                if(p<q){
                    Exp(pauliString, sign * angle, qubits);
                }
                else{
                    Exp(pauliString, -1.0 * sign * angle, qubits);
                }

            }
        }
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Applies time-evolution by a cluster operator PQQR term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a cluster operator PQQR term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerClusterOperatorPQQRTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let angle = (0.25 * coeff[0]);
            let p = idxFermions[0];
            let q = idxFermions[1];
            let r = idxFermions[3];

            // For all cases, do the same thing:
            // p < r < q (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p Y_r - Y_p X_r) (same as Hermitian conjugate of r < p < q)
            // q < p < r (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p Y_r - Y_p X_r)
            // p < q < r (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p Y_r - Y_p X_r)
            
            // This amounts to applying a PQ term, followed by same PQ term after a CNOT from q to the parity bit.
            let termPR0 = GeneratorIndex((idxTermType, [angle * 2.0]), [p, r]);
            _ApplyJordanWignerClusterOperatorPQTerm(termPR0, stepSize, qubits);

            let pNew = q > p ? p | p - 1;
            let rNew = q > r ? r | r - 1;
            let termPR1 = GeneratorIndex((idxTermType, [-1.0 * angle * 2.0]), [pNew, rNew]);
            _ApplyJordanWignerClusterOperatorPQTerm(termPR1, stepSize, Exclude([q], qubits));
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Applies time-evolution by a cluster operator PQRS term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a cluster operator PQRS term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerClusterOperatorPQRSTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let p = idxFermions[0];
            let q = idxFermions[1];
            let r = idxFermions[2];
            let s = idxFermions[3];
            let angle = 0.125 * coeff[0] * stepSize;
            
            let x = PauliX;
            let y = PauliY;

            let ops = [[y,y,x,y],[x,x,x,y],[x,y,y,y],[y,x,y,y],[x,y,x,x],[y,x,x,x],[y,y,y,x],[x,x,y,x]];
            let (sortedIndices, signs, globalsign) = _JordanWignerClusterOperatorPQRSTermSigns([p,q,r,s]);

            for (idxOp in 0 .. Length(ops) - 1) {
                let pauliString = _ComputeJordanWignerPauliString(Length(qubits), idxFermions, ops[idxOp]);
                let sign = signs[idxOp];
                Exp(pauliString, globalsign*sign * angle, qubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    function _JordanWignerClusterOperatorPQRSTermSigns(indices: Int[]) : (Int[], Double[], Double)
    {
        let p = indices[0];
        let q = indices[1];
        let r = indices[2];
        let s = indices[3];
        mutable sorted = new Int[4];
        mutable signs = new Double[8];
        mutable sign = 1.0;

        if(p>q){
            set sign = sign * -1.0;
        }
        if(r>s){
            set sign = sign * -1.0;
        }
        if(Min([p,q]) > Min([r,s])){
            set sign = sign * -1.0;
            set sorted = [Min([r,s]), Max([r,s]), Min([p,q]), Max([p,q])];
        }
        else{
            set sorted = [Min([p,q]), Max([p,q]), Min([r,s]), Max([r,s])];
        }
        // sorted is now in the order
        // [p`,q`,r`,s`], where p`<q`; r`<s`, and Min(p`,q`) is smaller than Min(r`,s`).

        let p1 = sorted[0];
        let q1 = sorted[1];
        let r1 = sorted[2];
        let s1 = sorted[3];

        // Case (p,q) < (r,s) and (p,q) > (r,s)
        if(q1 < r1){
            // p1 < q1 < r1 < s1
            return ([p1,q1,r1,s1],[1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0], sign);
        }
        // Case interleaved
        elif(q1 > r1 && q1 < s1){
            // p1 < r1 < q1 < s1
            return ([p1,r1,q1,s1],[-1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0], sign);
        }
        // Case contained
        elif(q1 > r1 && q1 > s1){
            // p1 < r1 < s1 < q1
            return ([p1,r1,s1,q1],[1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0], sign);
        }
        else{
            fail("Completely invalid cluster operator specified.");
        }
    }

    /// # Summary
    /// Converts a Hamiltonian described by `JWOptimizedHTerms`
    /// to a `GeneratorSystem` expressed in terms of the
    /// `GeneratorIndex` convention defined in this file.
    ///
    /// # Input
    /// ## data
    /// Description of Hamiltonian in `JWOptimizedHTerms` format.
    ///
    /// # Output
    /// Representation of Hamiltonian as `GeneratorSystem`.
    function JordanWignerClusterOperatorGeneratorSystem (data : JordanWignerInputState[]) : GeneratorSystem {
        return GeneratorSystem(Length(data), _JordanWignerClusterOperatorGeneratorSystemImpl(data, _));
    }

    function _JordanWignerClusterOperatorGeneratorSystemImpl(data : JordanWignerInputState[], idx: Int) : GeneratorIndex {
        return _JordanWignerClusterOperatorGeneratorIndex(data[idx]);
    }

    function _JordanWignerClusterOperatorGeneratorIndex(data: JordanWignerInputState): GeneratorIndex {
        let ((real, imaginary), idxFermions) = data!;
        if(Length(idxFermions) == 2){
            // PQ term
            return GeneratorIndex(([0],[real]),idxFermions);
        }
        elif(Length(idxFermions) == 4){
            if(idxFermions[1] == idxFermions[2]){
                // PQQR term
                return GeneratorIndex(([1],[real]),idxFermions);
            }
            else{
                // PQRS term
                return GeneratorIndex(([2],[real]),idxFermions);
            }
        }
        else{
            // Any other term in invalid
            return GeneratorIndex(([-1],[0.0]),[0]);
        }
    }

    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// See [Dynamical Generator Modeling](../libraries/data-structures#dynamical-generator-modeling)
    /// for more details.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the JordanWigner.
    /// ## stepSize
    /// Dummy variable to match signature of simulation algorithms.
    /// ## qubits
    /// Register acted upon by time-evolution operator.
    operation _JordanWignerClusterOperatorImpl(generatorIndex : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, idxDoubles), idxFermions) = generatorIndex!;
            let termType = idxTermType[0];
            
            if (termType == 0) {
                _ApplyJordanWignerClusterOperatorPQTerm (generatorIndex, stepSize, qubits);
            }
            elif (termType == 1) {
                _ApplyJordanWignerClusterOperatorPQQRTerm(generatorIndex, stepSize, qubits);
            }
            elif (termType == 2) {
                _ApplyJordanWignerClusterOperatorPQRSTerm (generatorIndex, stepSize, qubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the JordanWigner.
    ///
    /// # Output
    /// An `EvolutionUnitary` representing time-evolution by the term
    /// referenced in `generatorIndex.
    function _JordanWignerClusterOperatorFunction (generatorIndex : GeneratorIndex) : EvolutionUnitary {
        
        return EvolutionUnitary(_JordanWignerClusterOperatorImpl(generatorIndex, _, _));
    }
    
    
    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// # Output
    /// An `EvolutionSet` that maps a `GeneratorIndex` for the JordanWigner basis to
    /// an `EvolutionUnitary.
    function JordanWignerClusterOperatorEvolutionSet () : EvolutionSet {
        
        return EvolutionSet(_JordanWignerClusterOperatorFunction(_));
    }
    
}


