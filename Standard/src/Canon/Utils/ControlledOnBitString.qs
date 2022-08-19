// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Applies a unitary operation on the target register, controlled on a a state specified by a given
    /// bit mask.
    ///
    /// # Input
    /// ## bits
    /// The bit string to control the given unitary operation on.
    /// ## oracle
    /// The unitary operation to be applied on the target register.
    /// ## targetRegister
    /// The target register to be passed to `oracle` as an input.
    /// ## controlRegister
    /// A quantum register that controls application of `oracle`.
    ///
    /// # Remarks
    /// The pattern given by `bits` may be shorter than `controlRegister`,
    /// in which case additional control qubits are ignored (that is, neither
    /// controlled on $\ket{0}$ nor $\ket{1}$).
    /// If `bits` is longer than `controlRegister`, an error is raised.
    ///
    /// For example, `bits = [0,1,0,0,1]` means that `oracle` is applied if and only if `controlRegister`
    /// is in the state $\ket{0}\ket{1}\ket{0}\ket{0}\ket{1}$.
    operation ApplyControlledOnBitString<'T> (bits : Bool[], oracle : ('T => Unit is Adj + Ctl), controlRegister : Qubit[], targetRegister : 'T)
    : Unit is Adj + Ctl {
        // The control register must have enough bits to implement the requested control.
        Fact(Length(bits) <= Length(controlRegister), "Control register shorter than control pattern.");

        // Use a subregister of the controlled register when
        // bits is shorter than controlRegister.
        let controlSubregister = controlRegister[...Length(bits) - 1];
        within {
            ApplyPauliFromBitString(PauliX, false, bits, controlSubregister);
        } apply {
            Controlled oracle(controlSubregister, targetRegister);
        }
    }


    /// # Summary
    /// Returns a unitary operation that applies an oracle on the target register if the control register state corresponds to a specified bit mask.
    ///
    /// # Description
    /// The output of this function is an operation that can be represented by a
    /// unitary transformation $U$ such that
    /// \begin{align}
    ///     U \ket{b_0 b_1 \cdots b_{n - 1}} \ket{\psi} = \ket{b_0 b_1 \cdots b_{n-1}} \otimes
    ///     \begin{cases}
    ///         V \ket{\psi} & \textrm{if} (b_0 b_1 \cdots b_{n - 1}) = \texttt{bits} \\\\
    ///         \ket{\psi} & \textrm{otherwise}
    ///     \end{cases},
    /// \end{align}
    /// where $V$ is a unitary transformation that represents the action of the
    /// `oracle` operation.
    ///
    /// # Input
    /// ## bits
    /// The bit string to control the given unitary operation on.
    /// ## oracle
    /// The unitary operation to be applied on the target register.
    ///
    /// # Output
    /// A unitary operation that applies `oracle` on the target register if the control register state corresponds to the bit mask `bits`.
    ///
    /// # Remarks
    /// The pattern given by `bits` may be shorter than `controlRegister`,
    /// in which case additional control qubits are ignored (that is, neither
    /// controlled on $\ket{0}$ nor $\ket{1}$).
    /// If `bits` is longer than `controlRegister`, an error is raised.
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
    function ControlledOnBitString<'T> (bits : Bool[], oracle : ('T => Unit is Adj + Ctl))
    : ((Qubit[], 'T) => Unit is Adj + Ctl) {
        return ApplyControlledOnBitString(bits, oracle, _, _);
    }


    /// # Summary
    /// Applies a unitary operation on the target register if the control
    /// register state corresponds to a specified nonnegative integer.
    ///
    /// # Input
    /// ## numberState
    /// A nonnegative integer on which the operation `oracle` should be
    /// controlled.
    /// ## oracle
    /// A unitary operation to be controlled.
    /// ## targetRegister
    /// A register on which to apply `oracle`.
    /// ## controlRegister
    /// A quantum register that controls application of `oracle`.
    ///
    /// # Remarks
    /// The value of `numberState` is interpreted using a little-endian encoding.
    ///
    /// `numberState` must be at most $2^\texttt{Length(controlRegister)} - 1$.
    /// For example, `numberState = 537` means that `oracle`
    /// is applied if and only if `controlRegister` is in the state $\ket{537}$.
    operation ApplyControlledOnInt<'T> (numberState : Int, oracle : ('T => Unit is Adj + Ctl), controlRegister : Qubit[], targetRegister : 'T)
    : Unit is Adj + Ctl {
        let bits = IntAsBoolArray(numberState, Length(controlRegister));
        (ControlledOnBitString(bits, oracle))(controlRegister, targetRegister);
    }


    /// # Summary
    /// Returns a unitary operator that applies an oracle on the target register
    /// if the control register state corresponds to a specified nonnegative integer.
    ///
    /// # Input
    /// ## numberState
    /// Nonnegative integer.
    /// ## oracle
    /// Unitary operator.
    ///
    /// # Output
    /// A unitary operator that applies `oracle` on the target register if the
    /// control register state corresponds to the number state `numberState`.
    ///
    /// # Remarks
    /// The value of `numberState` is interpreted using a little-endian encoding.
    function ControlledOnInt<'T>(numberState : Int, oracle : ('T => Unit is Adj + Ctl))
    : ((Qubit[], 'T) => Unit is Adj + Ctl) {
        return ApplyControlledOnInt(numberState, oracle, _, _);
    }

}

