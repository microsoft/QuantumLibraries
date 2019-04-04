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
    /// Number of fermions in the system.
	/// ## fermionIndices
	/// fermionic operator indices.
    ///
    /// # Output
    /// Bitstring `Bool[]` that is `true` where a `PauliZ` should be applied.
	function _ComputeJordanWignerBitString_(nFermions: Int,  fermionIndices: Int[]) : Bool[] {
		if(Length(fermionIndices) % 2 != 0){
			fail $"ComputeJordanWignerString failed. `fermionIndices` must contain an even number of terms.";
		}
		
		mutable zString = new Bool[nFermions];
		for(fermionIdx in fermionIndices){
			if(fermionIdx >= nFermions){
				fail $"ComputeJordanWignerString failed. fermionIdx {fermionIdx} out of range.";
			}
			for(idx in 0..fermionIdx){
				set zString[idx] = zString[idx] ? false | true;
			}
		}
		
		for(fermionIdx in fermionIndices){
			set zString[fermionIdx] = false;
		}
		return zString;
	}

	// Identical to `_ComputeJordanWignerBitString_`, except with the map  
	// false -> PauliI and true -> PauliZ
	function _ComputeJordanWignerPauliZString_(nFermions: Int,  fermionIndices: Int[]) : Pauli[] {
		let bitString = _ComputeJordanWignerBitString_(nFermions, fermionIndices);

		mutable pauliString = new Pauli[nFermions];
		for(idx in 0..nFermions-1){
			set pauliString[idx] = bitString[idx] ? PauliZ | PauliI;
		}

		return pauliString;
	}

	// Identical to `_ComputeJordanWignerPauliZString_`, except that some
	// specified elements are substituted.
	function _ComputeJordanWignerPauliString_(nFermions: Int,  fermionIndices: Int[], pauliReplacements : Pauli[]) : Pauli[] {
		mutable pauliString = _ComputeJordanWignerPauliZString_(nFermions, fermionIndices);

		for(idx in 0..Length(fermionIndices)-1){
			let idxFermion = fermionIndices[idx];
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
	operation _ApplyJordanWignerClusterOperatorPQTerm_ (term : GeneratorIndex, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let p = idxFermions[0];
			let q = idxFermions[1];
			let angle = 0.5 * coeff[0];
			let ops = [[PauliX, PauliY], [PauliY, PauliX]];
            let signs = [+1.0, -1.0];

            for (idxOp in 0 .. Length(ops) - 1) {
				let pauliString = _ComputeJordanWignerPauliString_(Length(qubits), idxFermions, ops[idxOp]);
				let sign = signs[idxOp];
				if(p < q){
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
    operation _ApplyJordanWignerClusterOperatorPQQRTerm_ (term : GeneratorIndex, qubits : Qubit[]) : Unit {
        
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
            _ApplyJordanWignerClusterOperatorPQTerm_(termPR0, qubits);

			let pNew = q > p ? p | p - 1;
			let rNew = q > r ? r | r - 1;
			let termPR1 = GeneratorIndex((idxTermType, [-1.0 * angle * 2.0]), [pNew, rNew]);
			_ApplyJordanWignerClusterOperatorPQTerm_(termPR1, Exclude([q], qubits));
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
    /// `GeneratorIndex` representing a cluster operator PQQR term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerClusterOperatorPQRSTerm_ (term : GeneratorIndex, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let p = idxFermions[0];
			let q = idxFermions[1];
			let r = idxFermions[2];
			let s = idxFermions[3];
			let angle = 0.125 * coeff[0];
			
			let ops = [
			[PauliX, PauliX, PauliX, PauliX], 
			[PauliX, PauliX, PauliY, PauliY], 
			[PauliX, PauliY, PauliX, PauliY], 
			[PauliY, PauliX, PauliX, PauliY], 
			[PauliY, PauliY, PauliY, PauliY], 
			[PauliY, PauliY, PauliX, PauliX], 
			[PauliY, PauliX, PauliY, PauliX], 
			[PauliX, PauliY, PauliY, PauliX]];
            
			
			let ops = [[PauliX, PauliY], [PauliY, PauliX]];
            let signs = [+1.0, -1.0];

            for (idxOp in 0 .. Length(ops) - 1) {
				let pauliString = _ComputeJordanWignerPauliString_(Length(qubits), idxFermions, ops[idxOp]);
				let sign = signs[idxOp];
				if(p < q){
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
    
}


