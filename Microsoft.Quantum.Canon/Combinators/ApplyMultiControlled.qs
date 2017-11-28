namespace Microsoft.Quantum.Canon
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Combinators for constructing multiply controlled versions of operations
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /// # Summary 
    /// The signature type of CCNOT gate 
    newtype CCNOTop = (( Qubit, Qubit, Qubit ) => () : Adjoint );

    /// # Summary 
    /// Applies a multiply controlled version of singly controlled 
    /// operation.
    /// # Inputs
    /// ## singlyControlledOp
    /// Singly controlled operation. The first qubit in the argument of the operation 
    /// assumed to be a control and the rest are assumed to be target qubits.
    /// ## ccnot 
    /// The ccnot gate to use for the construction
    /// ApplyMultiControlled calls singlyControlledOp with the argument of the length at least 1.
    /// ## controls
    /// The qubits an operation being controlled on, the length of controls must be at least 1.
    /// ## targets
    /// The target qubits
    /// # Remarks 
    /// This operation uses only clean ancilla qubits. 
    /// # References
    /// - [ *Michael A. Nielsen , Isaac L. Chuang*,
    ///     Quantum Computation and Quantum Information ](http://doi.org/10.1017/CBO9780511976667)
    /// # See Also 
    /// - For the explanation and circuit diagram see Figure 4.10, Section 4.3 in Nielsen & Chuang
    operation ApplyMultiControlledC (
        singlyControlledOp : ( Qubit[] => () ),
        ccnot : CCNOTop,
        controls : Qubit[],
        targets : Qubit[] ) : () {
        body {
            AssertBoolEqual(
                Length(controls) >= 1, true,
                "Length of controls must be at least 1" );
            if( Length(controls) == 1 ) {
                singlyControlledOp( controls + targets );
            } else {
                using( ancillas = Qubit[ Length(controls) - 1 ] ) {
                    AndLadder(ccnot, controls, ancillas);
                    singlyControlledOp( [Tail(ancillas)] + targets );
                    (Adjoint AndLadder)(ccnot, controls, ancillas);
                }
            }
        }
        controlled( extraControls ) {
            ApplyMultiControlledC(singlyControlledOp,ccnot,extraControls+controls,targets);
        }
    }

    /// # See Also
    /// - @"Microsoft.Quantum.Canon.ApplyMultiControlledCA"
    operation ApplyMultiControlledCA (
        singlyControlledOp : ( Qubit[] => () : Adjoint ),
        ccnot : CCNOTop,
        controls : Qubit[],
        targets : Qubit[] ) : () {
        body {
            AssertBoolEqual(
                Length(controls) >= 1, true,
                "Length of controls must be at least 1" );
            if( Length(controls) == 1 ) {
                singlyControlledOp( controls + targets );
            } else {
                using( ancillas = Qubit[ Length(controls) - 1 ] ) {
                    AndLadder(ccnot, controls, ancillas);
                    singlyControlledOp( [Tail(ancillas)] + targets );
                    (Adjoint AndLadder)(ccnot, controls, ancillas);
                }
            }
        }
        adjoint auto
        controlled( extraControls ) {
            ApplyMultiControlledCA(singlyControlledOp,ccnot,extraControls+controls,targets);
        }
        controlled adjoint auto
    }

    /// # Inputs 
    /// ## ccnot
    /// The CCNOT gate to use for the construction
    /// ## controls
    /// We use notation |x₁,…,xₙ⟩ for values of these qubit in computational basis.
    /// These qubits are left unchanged. 
    /// The length of controls must be at least 2 and equal to one plus the length of targets .
    /// ## targets 
    /// We use notation |x₁,…,xₙ⟩ for values of these qubit in computational basis.
    /// The length of targets must be at least 1 and equal to the length of controls minus one. 
    /// # Summary 
    /// Applies a unitary given by the following map on computational basis vectors:
    /// |x₁,…,xₙ⟩|y₁,…,yₙ₋₁⟩ ↦ |x₁,…,xₙ⟩|y₁⊕(x₁∧x₂),…,yₙ₋₁⊕(x₁∧x₂∧…∧xₙ)⟩
    /// # References
    /// - [ *Michael A. Nielsen , Isaac L. Chuang*,
    ///     Quantum Computation and Quantum Information ](http://doi.org/10.1017/CBO9780511976667)
    /// # See Also 
    /// - Used as a part of @"Microsoft.Quantum.Canon.ApplyMultiControlledC"
    ///   and  @"Microsoft.Quantum.Canon.ApplyMultiControlledCA"
    /// - For the explanation and circuit diagram see Figure 4.10, Section 4.3 in Nielsen & Chuang
    operation AndLadder ( ccnot : CCNOTop, controls : Qubit[], targets  : Qubit[]) : () {
		body{
            AssertBoolEqual( Length(controls) == Length(targets) + 1, true,
                "Length(controls) must be equal to Length(target) + 1" );
            AssertBoolEqual( Length(controls) >= 2, true,
                "The operation is not defined for less than 2 controls" );
			ccnot(controls[0], controls[1], targets[0]);
			for ( k in 1 .. Length(targets)- 1 ) {
				ccnot(controls[k + 1],targets[k - 1],targets[k]);
			}
		}
		adjoint auto
	}
}
