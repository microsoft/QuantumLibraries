// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Arrays;
    
    
    // This block encoding for qubitization runs off data optimized for a jordan-wigner encoding.
    // This collects terms Z, ZZ, PQandPQQR, hpqrs separately.
    // This only apples the needed hpqrs XXXX XXYY terms.
    
    // Convention for GeneratorIndex = ((Int[],Double[]), Int[])
    // We index single Paulis as 0 for I, 1 for X, 2 for Y, 3 for Z.
    // We index Pauli strings with arrays of integers e.g. a = [3,1,1,2] for ZXXY.
    // We assume the zeroth element of Double[] is the angle of rotation
    // We index the qubits that Pauli strings act on with arrays of integers e.g. q = [2,4,5,8] for Z_2 X_4 X_5, Y_8
    // An example of a Pauli string GeneratorIndex is thus ((a,b), q)
    
    // Consider the Hamiltonian H = 0.1 XI + 0.2 IX + 0.3 ZY
    // Its GeneratorTerms are (([1],b),[0]), 0.1),  (([1],b),[1]), 0.2),  (([3,2],b),[0,1]), 0.3).
    
    /// # Summary
    /// Converts a GeneratorIndex describing a Z term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a Z term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing Z term as Pauli terms.
    function _ZTermToPauliGenIdx(term : GeneratorIndex) : GeneratorIndex[] {
        let ((idxTermType, coeff), idxFermions) = term!;
        return [GeneratorIndex(([3], coeff), idxFermions)];
    }
    
    
    /// # Summary
    /// Converts a GeneratorIndex describing a ZZ term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a ZZ term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing ZZ term as Pauli terms.
    function _ZZTermToPauliGenIdx (term : GeneratorIndex) : GeneratorIndex[] {
        let ((idxTermType, coeff), idxFermions) = term!;
        return [GeneratorIndex(([3, 3], coeff), idxFermions)];
    }
    
    
    /// # Summary
    /// Converts a GeneratorIndex describing a PQ term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQ term as Pauli terms.
    function _PQTermToPauliGenIdx_ (term : GeneratorIndex) : GeneratorIndex[] {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        let newCoeff = [coeff[0]];
        let qubitPidx = idxFermions[0];
        let qubitQidx = idxFermions[1];
        let qubitIndices = IntArrayFromRange(qubitPidx .. qubitQidx);
        return [GeneratorIndex((([1] + ConstantArray(Length(qubitIndices) - 2, 3)) + [1], newCoeff), qubitIndices), GeneratorIndex((([2] + ConstantArray(Length(qubitIndices) - 2, 3)) + [2], newCoeff), qubitIndices)];
    }
    
    
    /// # Summary
    /// Converts a GeneratorIndex describing a PQ or PQQR term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ or PQQR term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQ or PQQR term as Pauli terms.
    function _PQandPQQRTermToPauliGenIdx_ (term : GeneratorIndex) : GeneratorIndex[] {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        let newCoeff = [coeff[0]];
        
        if (Length(idxFermions) == 2) {
            return _PQTermToPauliGenIdx_(term);
        }
        else {
            let qubitPidx = idxFermions[0];
            let qubitQidx = idxFermions[1];
            let qubitRidx = idxFermions[3];
            
            if (qubitPidx < qubitQidx and qubitQidx < qubitRidx) {
                
                // Apply XZ..ZIZ..ZX
                let qubitIndices = IntArrayFromRange(qubitPidx .. qubitQidx - 1) + IntArrayFromRange(qubitQidx + 1 .. qubitRidx);
                return [GeneratorIndex((([1] + ConstantArray(Length(qubitIndices) - 2, 3)) + [1], newCoeff), qubitIndices), GeneratorIndex((([2] + ConstantArray(Length(qubitIndices) - 2, 3)) + [2], newCoeff), qubitIndices)];
            }
            else {
                
                // Apply ZI..IXZ..ZX or XZ..ZXI..IZ
                let qubitIndices = IntArrayFromRange(qubitPidx .. qubitRidx) + [qubitQidx];
                return [GeneratorIndex((([1] + ConstantArray(Length(qubitIndices) - 3, 3)) + [1, 3], newCoeff), qubitIndices), GeneratorIndex((([2] + ConstantArray(Length(qubitIndices) - 3, 3)) + [2, 3], newCoeff), qubitIndices)];
            }
        }
    }
    
    
    /// # Summary
    /// Converts a GeneratorIndex describing a PQRS term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQRS term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQRS term as Pauli terms.
    function _V0123TermToPauliGenIdx_ (term : GeneratorIndex) : GeneratorIndex[] {
        
        let ((idxTermType, v0123), idxFermions) = term!;
        let qubitsPQ = idxFermions[0 .. 1];
        let qubitsRS = idxFermions[2 .. 3];
        let qubitsPQJW = IntArrayFromRange(qubitsPQ[0] + 1 .. qubitsPQ[1] - 1);
        let qubitsRSJW = IntArrayFromRange(qubitsRS[0] + 1 .. qubitsRS[1] - 1);
        let ops = [[1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], [2, 1, 1, 2], [2, 2, 2, 2], [2, 2, 1, 1], [2, 1, 2, 1], [1, 2, 2, 1]];
        mutable genIdxes = new GeneratorIndex[8];
        mutable nonZero = 0;
        
        for (idxOp in IndexRange(ops)) {
            
            if (IsNotZero(v0123[idxOp % 4])) {
                let newCoeff = [v0123[idxOp % 4]];
                set genIdxes[nonZero] = GeneratorIndex((ops[idxOp] + ConstantArray(Length(qubitsPQJW) + Length(qubitsRSJW), 3), newCoeff), ((qubitsPQ + qubitsRS) + qubitsPQJW) + qubitsRSJW);
                set nonZero = nonZero + 1;
            }
        }
        
        return genIdxes[0 .. nonZero - 1];
    }
    
    
    /// # Summary
    /// Converts a Hamiltonian described by `JWOptimizedHTerms`
    /// to a `GeneratorSystem` expressed in terms of the Pauli
    /// `GeneratorIndex`.
    ///
    /// # Input
    /// ## data
    /// Description of Hamiltonian in `JWOptimizedHTerms` format.
    ///
    /// # Output
    /// Representation of Hamiltonian as `GeneratorSystem`.
    function JordanWignerBlockEncodingGeneratorSystem (data : JWOptimizedHTerms) : GeneratorSystem {
        
        let (ZData, ZZData, PQandPQQRData, h0123Data) = data!;
        mutable genIdxes = new GeneratorIndex[((Length(ZData) + Length(ZZData)) + 2 * Length(PQandPQQRData)) + 8 * Length(h0123Data)];
        mutable startIdx = 0;
        
        for (idx in IndexRange(ZData)) {
            
            // Array of Arrays of Length 1
            set genIdxes[idx] = (_ZTermToPauliGenIdx(HTermToGenIdx(ZData[idx], [0])))[0];
        }
        
        set startIdx = Length(ZData);
        
        for (idx in IndexRange(ZZData)) {
            
            // Array of Arrays of Length 1
            set genIdxes[startIdx + idx] = (_ZZTermToPauliGenIdx(HTermToGenIdx(ZZData[idx], [1])))[0];
        }
        
        set startIdx = startIdx + Length(ZZData);
        
        for (idx in IndexRange(PQandPQQRData)) {
            
            // Array of Arrays of Length 2
            let genArr = _PQandPQQRTermToPauliGenIdx_(HTermToGenIdx(PQandPQQRData[idx], [2]));
            set genIdxes[startIdx + 2 * idx] = genArr[0];
            set genIdxes[(startIdx + 2 * idx) + 1] = genArr[1];
        }
        
        set startIdx = startIdx + 2 * Length(PQandPQQRData);
        mutable finalIdx = startIdx;
        
        for (idx in 0 .. Length(h0123Data) - 1) {
            
            // Array of Arrays of Length up to 8
            let genArr = _V0123TermToPauliGenIdx_(HTermToGenIdx(h0123Data[idx], [3]));
            
            for (idx0123 in IndexRange(genArr)) {
                set genIdxes[finalIdx] = genArr[idx0123];
                set finalIdx = finalIdx + 1;
            }
        }
        
        return GeneratorSystem(finalIdx, LookupFunction(genIdxes[0 .. finalIdx - 1]));
    }
    
}


