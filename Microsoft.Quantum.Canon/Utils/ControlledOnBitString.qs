// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Applies a unitary operator on the target register if the control register state corresponds to a specified bit mask.
    ///
    /// # Input
    /// ## bits
    /// Boolean array.
    /// ## oracle
    /// Unitary operator.
    /// ## targetRegister
    /// Quantum register acted on by `oracle`.
    /// ## controlRegister
    /// Quantum register that controls application of `oracle`.
    ///
    /// # Remarks 
    /// The length of `bits` and `controlRegister` must be equal.
    /// For example, `bits = [0;1;0;0;1]` means that `oracle` is applied if and only if `controlRegister`" is in the state $\ket{0}\ket{1}\ket{0}\ket{0}\ket{1}$.
    operation ControlledOnBitStringImpl(bits : Bool[] , oracle: (Qubit[] => (): Adjoint, Controlled), controlRegister : Qubit[], targetRegister: Qubit[]) : ()
    {
	    body{
		    WithCA(ApplyPauliFromBitString(PauliX, false, bits, _), (Controlled oracle)(_, targetRegister), controlRegister);
	    }

	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }

    /// # Summary
    /// Returns a unitary operator that applies an oracle on the target register if the control register state corresponds to a specified bit mask.
    ///
    /// # Input
    /// ## bits
    /// Boolean array.
    /// ## oracle
    /// Unitary operator.
    ///
    /// # Output
    /// A unitary operator that applies `oracle` on the target register if the control register state corresponds to the bit mask `bits`.
    ///
    /// # Remarks 
    /// Obtained by partial application of @"microsoft.quantum.canon.ControlledOnBitStringImpl".
    function ControlledOnBitString(bits : Bool[] , oracle: (Qubit[] => (): Adjoint, Controlled)) : ((Qubit[],Qubit[]) => (): Adjoint, Controlled)
    {
	    return ControlledOnBitStringImpl(bits, oracle, _, _);
    }

    /// # Summary
    /// Applies a unitary operator on the target register if the control register state corresponds to a specified positive integer.
    ///
    /// # Input
    /// ## numberState
    /// Positive integer.
    /// ## oracle
    /// Unitary operator.
    /// ## targetRegister
    /// Quantum register acted on by `oracle`.
    /// ## controlRegister
    /// Quantum register that controls application of `oracle`.
    ///
    /// # Remarks 
    /// `numberState` must be at most $2^\texttt{Length(controlRegister)} - 1$.
    /// For example, `numberState = 537` means that `oracle` is applied if and only if `controlRegister` is in the state $\ket{537}$.
    operation ControlledOnIntImpl(numberState : Int , oracle: (Qubit[] => (): Adjoint, Controlled), controlRegister : Qubit[], targetRegister: Qubit[]) : ()
    {
	    body {

		    let bits = BoolArrFromPositiveInt(numberState, Length(controlRegister));

		    (ControlledOnBitString(bits, oracle))(controlRegister, targetRegister);

	    }

	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }

    /// # Summary
    /// Returns a unitary operator that applies an oracle on the target register if the control register state corresponds to a specified positive integer.
    ///
    /// # Input
    /// ## numberState
    /// Positive integer.
    /// ## oracle
    /// Unitary operator.
    ///
    /// # Output
    /// A unitary operator that applies `oracle` on the target register if the control register state corresponds to the number state `numberState`.
    ///
    /// # Remarks 
    /// Obtained by partial application of @"microsoft.quantum.canon.ControlledOnIntImpl".
    function ControlledOnInt(numberState : Int , oracle: (Qubit[] => (): Adjoint, Controlled)) : ((Qubit[],Qubit[]) => (): Adjoint, Controlled)
    {
	    return ControlledOnIntImpl(numberState, oracle, _, _);
    }

}
