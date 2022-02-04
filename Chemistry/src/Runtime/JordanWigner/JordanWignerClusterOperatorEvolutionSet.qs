// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Computes Z component of Jordan–Wigner string between
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
        if Length(idxFermions) % 2 != 0 {
            fail $"ComputeJordanWignerString failed. `idxFermions` must contain an even number of terms.";
        }

        mutable zString = [false, size = nFermions];
        for fermionIdx in idxFermions {
            if fermionIdx >= nFermions {
                fail $"ComputeJordanWignerString failed. fermionIdx {fermionIdx} out of range.";
            }
            for idx in 0..fermionIdx {
                set zString w/= idx <- not zString[idx];
            }
        }

        for fermionIdx in idxFermions {
            set zString w/= fermionIdx <- false;
        }
        return zString;
    }

    // Identical to `_ComputeJordanWignerBitString`, except with the map  
    // false -> PauliI and true -> PauliZ
    function _ComputeJordanWignerPauliZString(nFermions: Int,  idxFermions: Int[]) : Pauli[] {
        let bitString = _ComputeJordanWignerBitString(nFermions, idxFermions);
        return BoolArrayAsPauli (PauliZ, true, bitString);
    }

    // Identical to `_ComputeJordanWignerPauliZString`, except that some
    // specified elements are substituted.
    function _ComputeJordanWignerPauliString(nFermions: Int,  idxFermions: Int[], pauliReplacements : Pauli[]) : Pauli[] {
        mutable pauliString = _ComputeJordanWignerPauliZString(nFermions, idxFermions);

        for idx in IndexRange(idxFermions) {
            let idxFermion = idxFermions[idx];
            let op = pauliReplacements[idx];
            set pauliString w/= idxFermion <- op;
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
    internal operation _ApplyJordanWignerClusterOperatorPQTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, coeff), idxFermions) = term!;
        let p = idxFermions[0];
        let q = idxFermions[1];
        if p == q {
            fail $"Unitary coupled-cluster PQ failed: indices {p}, {q} must be distinct";
        }
        let angle = 0.5 * coeff[0] * stepSize;
        let ops = [[PauliX, PauliY], [PauliY, PauliX]];
        let signs = [+1.0, -1.0];

        for (op, sign) in Zipped(ops, signs) {
            let pauliString = _ComputeJordanWignerPauliString(Length(qubits), idxFermions, op);
            Exp(pauliString, sign * angle, qubits);
        }
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
    internal operation _ApplyJordanWignerClusterOperatorPQRSTerm (term : GeneratorIndex, stepSize: Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, coeff), idxFermions) = term!;
        let p = idxFermions[0];
        let q = idxFermions[1];
        let r = idxFermions[2];
        let s = idxFermions[3];
        let angle = 0.125 * coeff[0] * stepSize;

        if p == q or p==r or p==s or q==r or q==s or r==s {
            fail($"Unitary coupled-cluster PQRS failed: indices {p}, {q}, {r}, {s} must be distinct");
        }

        let x = PauliX;
        let y = PauliY;

        let ops = [[y,y,x,y],[x,x,x,y],[x,y,y,y],[y,x,y,y],[x,y,x,x],[y,x,x,x],[y,y,y,x],[x,x,y,x]];
        let (sortedIndices, signs, globalSign) = _JordanWignerClusterOperatorPQRSTermSigns([p,q,r,s]);

        for (op, sign) in Zipped(ops, signs) {
            let pauliString = _ComputeJordanWignerPauliString(Length(qubits), sortedIndices, op);
            Exp(pauliString, globalSign * sign * angle, qubits);
        }
    }
    
    function _JordanWignerClusterOperatorPQRSTermSigns(indices: Int[]) : (Int[], Double[], Double) {
        let p = indices[0];
        let q = indices[1];
        let r = indices[2];
        let s = indices[3];
        mutable sorted = [0, size = 4];
        mutable signs = [0.0, size = 8];
        mutable sign = 1.0;

        if(p>q){
            set sign = sign * -1.0;
        }
        if(r>s){
            set sign = sign * -1.0;
        }
        if( Min([p,q]) > Min([r,s]) ){
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
        elif(q1 > r1 and q1 < s1){
            // p1 < r1 < q1 < s1
            return ([p1,r1,q1,s1],[-1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0], sign);
        }
        // Case contained
        elif(q1 > r1 and q1 > s1){
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
        return GeneratorSystem(Length(data), JordanWignerStateAsGeneratorIndex(data, _));
    }

    internal function JordanWignerStateAsGeneratorIndex(data : JordanWignerInputState[], idx : Int) : GeneratorIndex {
        let ((real, imaginary), idxFermions) = data[idx]!;

        if Length(idxFermions) == 2 {
            // PQ term
            return GeneratorIndex(([0], [real]), idxFermions);
        } elif Length(idxFermions) == 4 {
            // PQRS term
            return GeneratorIndex(([2], [real]), idxFermions);
        } else {
            // Any other term in invalid
            return GeneratorIndex(([-1], [0.0]), [0]);
        }
    }

    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// See [Dynamical Generator Modeling](xref:microsoft.quantum.libraries.overview.data-structures#dynamical-generator-modeling)
    /// for more details.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the JordanWigner.
    /// ## stepSize
    /// Dummy variable to match signature of simulation algorithms.
    /// ## qubits
    /// Register acted upon by time-evolution operator.
    internal operation _JordanWignerClusterOperatorImpl(generatorIndex : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let ((idxTermType, idxDoubles), idxFermions) = generatorIndex!;
        let termType = idxTermType[0];

        if termType == 0 {
            _ApplyJordanWignerClusterOperatorPQTerm (generatorIndex, stepSize, qubits);
        }
        elif termType == 2 {
            _ApplyJordanWignerClusterOperatorPQRSTerm (generatorIndex, stepSize, qubits);
        }
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
    internal function _JordanWignerClusterOperatorFunction(generatorIndex : GeneratorIndex) : EvolutionUnitary {
        return EvolutionUnitary(_JordanWignerClusterOperatorImpl(generatorIndex, _, _));
    }


    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// # Output
    /// An `EvolutionSet` that maps a `GeneratorIndex` for the JordanWigner basis to
    /// an `EvolutionUnitary.
    function JordanWignerClusterOperatorEvolutionSet() : EvolutionSet {
        return EvolutionSet(_JordanWignerClusterOperatorFunction);
    }
    
}

