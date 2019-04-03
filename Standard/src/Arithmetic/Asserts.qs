// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
	/// Asserts that the probability of a specific state of a quantum register has the
	/// expected value.
	///
    /// Given an $n$-qubit quantum state $\ket{\psi}=\sum^{2^n-1}_{j=0}\alpha_j \ket{j}$,
    /// asserts that the probability $|\alpha_j|^2$ of the state $\ket{j}$ indexed by $j$
    /// has the expected value.
    ///
    /// # Input
    /// ## stateIndex
    /// The index $j$ of the state $\ket{j}$ represented by a `LittleEndian`
    /// register.
    ///
    /// ## expected
    /// The expected value of $|\alpha_j|^2$.
    ///
    /// ## qubits
    /// The qubit register that stores $\ket{\psi}$ in little-endian format.
    ///
    /// ## tolerance
    /// Absolute tolerance on the difference between actual and expected.
    ///
    /// # Remarks
    /// ## Example
    /// Suppose that the `qubits` register encodes a 3-qubit quantum state
    /// $\ket{\psi}=\sqrt{1/8}\ket{0}+\sqrt{7/8}\ket{6}$ in little-endian format.
    /// This means that the number states $\ket{0}\equiv\ket{0}\ket{0}\ket{0}$
    /// and $\ket{6}\equiv\ket{0}\ket{1}\ket{1}$. Then the following asserts succeed:
    /// - `AssertProbInt(0,0.125,qubits,10e-10);`
    /// - `AssertProbInt(6,0.875,qubits,10e-10);`
    operation AssertProbIntLE(stateIndex : Int, expected : Double, qubits : LittleEndian, tolerance : Double) : Unit {
        let nQubits = Length(qubits!);
        let bits = BoolArrFromPositiveInt(stateIndex, nQubits);

        using (flag = Qubit()) {
            (ControlledOnBitString(bits, X))(qubits!, flag);
            AssertProb([PauliZ], [flag], One, expected, $"AssertProbInt failed on stateIndex {stateIndex}, expected probability {expected}.", tolerance);

            // Uncompute flag qubit.
            (ControlledOnBitString(bits, X))(qubits!, flag);
            Reset(flag);
        }
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertProbIntLE".
    operation AssertProbInt(stateIndex : Int, expected : Double, qubits : LittleEndian, tolerance : Double) : Unit {
        Renamed("Microsoft.Quantum.Arithmetic.AssertProbInt", "Microsoft.Quantum.Arithmetic.AssertProbIntLE");
        AssertProbIntLE(stateIndex, expected, qubits, tolerance);
    }

    /// # Summary
	/// Asserts that the probability of a specific state of a quantum register has the
	/// expected value.
    ///
    /// Given an $n$-qubit quantum state $\ket{\psi}=\sum^{2^n-1}_{j=0}\alpha_j \ket{j}$,
    /// asserts that the probability $|\alpha_j|^2$ of the state $\ket{j}$ indexed by $j$
    /// has the expected value.
    ///
    /// # Input
    /// ## stateIndex
    /// The index $j$ of the state $\ket{j}$ represented by a `BigEndian`
    /// register.
    ///
    /// ## expected
    /// The expected value of $|\alpha_j|^2$.
    ///
    /// ## qubits
    /// The qubit register that stores $\ket{\psi}$ in big-endian format.
    ///
    /// ## tolerance
    /// Absolute tolerance on the difference between actual and expected.
    ///
    /// # Remarks
    /// ## Example
    /// Suppose that the `qubits` register encodes a 3-qubit quantum state
    /// $\ket{\psi}=\sqrt{1/8}\ket{0}+\sqrt{7/8}\ket{6}$ in big-endian format.
    /// This means that the number states $\ket{0}\equiv\ket{0}\ket{0}\ket{0}$
    /// and $\ket{6}\equiv\ket{1}\ket{1}\ket{0}$. Then the following asserts succeed:
    /// - `AssertProbIntBE(0,0.125,qubits,10e-10);`
    /// - `AssertProbIntBE(6,0.875,qubits,10e-10);`
    operation AssertProbIntBE (stateIndex : Int, prob : Double, qubits : BigEndian, tolerance : Double) : Unit {
        AssertProbInt(stateIndex, prob, BigEndianAsLittleEndian(qubits), tolerance);
    }

    /// # Summary
    /// Asserts that the most significant qubit of a qubit register
    /// representing an unsigned integer is in a particular state.
    ///
    /// # Input
    /// ## value
    /// The value of the highest bit being asserted.
    /// ## number
    /// Unsigned integer of which the highest bit is checked.
    ///
    /// # Remarks
    /// The controlled version of this operation ignores controls.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Primitive.Assert
    operation AssertMostSignificantBitLE(value : Result, number : LittleEndian) : Unit {
        body (...) {
            let mostSingificantQubit = Tail(number!);
            Assert([PauliZ], [mostSingificantQubit], value, $"Most significant bit expected to be {value}");
        }

        adjoint self;

        controlled (ctrls, ...) {
            AssertHighestBit(value, number);
        }

        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertMostSignificantBitLE".
    operation AssertHighestBit(value : Result, number : LittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Arithmetic.AssertHighestBit", "Microsoft.Quantum.Arithmetic.AssertMostSignificantBitLE");
            AssertMostSignificantBitLE(value, number);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Summary
    /// Asserts that the `number` encoded in PhaseLittleEndian is less than `value`.
    ///
    /// # Input
    /// ## value
    /// `number` must be less than this.
    /// ## number
    /// Unsigned integer which is compared to `value`.
    ///
    /// # Remarks
    /// The controlled version of this operation ignores controls.
    operation AssertLessThanPhaseLE (value : Int, number : PhaseLittleEndian) : Unit
    {
        body (...)
        {
            let inner = ApplyLEOperationOnPhaseLEA(AssertHighestBit(One, _), _);
            WithA(Adjoint IntegerIncrementPhaseLE(value, _), inner, number);
        }
        
        adjoint self;
        
        controlled (ctrls, ...)
        {
            AssertLessThanPhaseLE(value, number);
        }
        
        controlled adjoint auto;
    }
    

}
