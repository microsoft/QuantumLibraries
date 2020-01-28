// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    
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
    /// For example, `bits = [0,1,0,0,1]` means that `oracle` is applied if and only if `controlRegister`" is in the state $\ket{0}\ket{1}\ket{0}\ket{0}\ket{1}$.
    operation ControlledOnBitStringImpl<'T> (bits : Bool[], oracle : ('T => Unit is Adj + Ctl), controlRegister : Qubit[], targetRegister : 'T) : Unit
    {
        body (...)
        {
            ApplyWithCA(ApplyPauliFromBitString(PauliX, false, bits, _), Controlled oracle(_, targetRegister), controlRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Returns a unitary operator that applies an oracle on the target register if the control register state corresponds to a specified bit mask.
    ///
    /// # Description
    /// The output of this function can be represented by a unitary transformation $U$ such that
    /// \begin{align}
    ///     U \ket{b_0 b_1 \cdots b_{n - 1}} \ket{\psi} = \ket{b_0 b_1 \cdots b_{n-1}} \otimes 
    ///     \begin{cases}
    ///         V \ket{\psi} & \textrm{if} (b_0 b_1 \cdots b_{n - 1}) = \texttt{bits} \\\\
    ///         \ket{\psi} & \textrm{otherwise}
    ///     \end{cases},
    /// \end{align}
    /// where $V$ is a unitary transformation that represents the action of the `oracle` operation.
    ///
    /// # Input
    /// ## bits
    /// The bit string to control the given unitary operator on.
    /// ## oracle
    /// Unitary operator to be applied on the target register.
    ///
    /// # Output
    /// A unitary operator that applies `oracle` on the target register if the control register state corresponds to the bit mask `bits`.
    ///
    /// # Remarks
    /// The length of `bits` and `controlRegister` must be equal.
    /// 
    /// Given a Boolean array `bits` and a unitary operation `oracle`, the output of this function
    /// is an operation that performs the following steps:
    /// * apply an `X` operation to each qubit of the control register that corresponds to `false` element of the `bits`;
    /// * apply `Controlled oracle` to the control and target registers;
    /// * apply an `X` operation to each qubit of the control register that corresponds to `false` element of the `bits` again to return the control register to the original state.
    ///
    /// The output of the `Controlled` functor is a special case of `ControlledOnBitString` where `bits` is equal to `[true, ..., true]`.
    ///
    /// # Example
    /// The following code snippets are equivalent:
    /// ```qsharp
    /// (ControlledOnBitString(bits, oracle))(controlRegister, targetRegister);
    /// ```
    /// and
    /// ```qsharp
    /// within {
    ///     ApplyPauliFromBitString(PauliX, false, bits, controlRegister);
    /// } apply {
    ///     Controlled oracle(controlRegister, targetRegister);
    /// }
    /// ```
    ///
    /// The following code prepares a state $\frac{1}{2}(\ket{00} - \ket{01} + \ket{10} + \ket{11})$:
    /// ```qsharp
    /// using (register = Qubit[2]) {
    ///     ApplyToEach(H, register);
    ///     (ControlledOnBitString([false], Z))(register[0..0], register[1]);
    /// }
    /// ```
    function ControlledOnBitString<'T> (bits : Bool[], oracle : ('T => Unit is Adj + Ctl)) : ((Qubit[], 'T) => Unit is Adj + Ctl)
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
    operation ControlledOnIntImpl<'T> (numberState : Int, oracle : ('T => Unit is Adj + Ctl), controlRegister : Qubit[], targetRegister : 'T) : Unit
    {
        body (...)
        {
            let bits = IntAsBoolArray(numberState, Length(controlRegister));
            (ControlledOnBitString(bits, oracle))(controlRegister, targetRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    function ControlledOnInt<'T> (numberState : Int, oracle : ('T => Unit is Adj + Ctl)) : ((Qubit[], 'T) => Unit is Adj + Ctl)
    {
        return ControlledOnIntImpl(numberState, oracle, _, _);
    }
    
}


