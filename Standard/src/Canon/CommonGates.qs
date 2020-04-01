// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies the controlled-X (CX) gate to a pair of qubits.
    ///
    /// $$
    /// \begin{align}
    ///     1 & 0 & 0 & 0 \\\\
    ///     0 & 1 & 0 & 0 \\\\
    ///     0 & 0 & 0 & 1 \\\\
    ///     0 & 0 & 1 & 0
    /// \end{align},
    /// $$
    /// where rows and columns are organized as in the quantum concepts guide.
    ///
    /// # Input
    /// ## control
    /// Control qubit for the CX gate.
    /// ## target
    /// Target qubit for the CX gate.
    ///
    /// # Remarks
    /// Equivalent to:
    /// ```qsharp
    /// Controlled X([control], target);
    /// ```
    /// and to:
    /// ```qsharp
    /// CNOT(control, target);
    /// ```
    operation CX(control : Qubit, target : Qubit) : Unit {
        body (...) {
            Controlled X([control], target);
        }
        adjoint self;
        controlled distribute;
        controlled adjoint self;
    }

    /// # Summary
    /// Applies the controlled-Y (CY) gate to a pair of qubits.
    ///
    /// $$
    /// \begin{align}
    ///     1 & 0 & 0 & 0 \\\\
    ///     0 & 1 & 0 & 0 \\\\
    ///     0 & 0 & 0 & -i \\\\
    ///     0 & 0 & i & 0
    /// \end{align},
    /// $$
    /// where rows and columns are organized as in the quantum concepts guide.
    ///
    /// # Input
    /// ## control
    /// Control qubit for the CY gate.
    /// ## target
    /// Target qubit for the CY gate.
    ///
    /// # Remarks
    /// Equivalent to:
    /// ```qsharp
    /// Controlled Y([control], target);
    /// ```
    operation CY(control : Qubit, target : Qubit) : Unit {
        body (...) {
            Controlled Y([control], target);
        }
        adjoint self;
        controlled distribute;
        controlled adjoint self;
    }

    /// # Summary
    /// Applies the controlled-Z (CZ) gate to a pair of qubits.
    ///
    /// $$
    /// \begin{align}
    ///     1 & 0 & 0 & 0 \\\\
    ///     0 & 1 & 0 & 0 \\\\
    ///     0 & 0 & 1 & 0 \\\\
    ///     0 & 0 & 0 & -1
    /// \end{align},
    /// $$
    /// where rows and columns are organized as in the quantum concepts guide.
    ///
    /// # Input
    /// ## control
    /// Control qubit for the CZ gate.
    /// ## target
    /// Target qubit for the CZ gate.
    ///
    /// # Remarks
    /// Equivalent to:
    /// ```qsharp
    /// Controlled Z([control], target);
    /// ```
    operation CZ(control : Qubit, target : Qubit) : Unit {
        body (...) {
            Controlled Z([control], target);
        }
        adjoint self;
        controlled distribute;
        controlled adjoint self;
    }

    /// # Summary
    /// Performs the Quantum Fourier Transform on a quantum register containing an
    /// integer in the big-endian representation.
    ///
    /// # Input
    /// ## qs
    /// Quantum register to which the Quantum Fourier Transform is applied
    ///
    /// # Remarks
    /// The input and output are assumed to be in big endian encoding.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyQuantumFourierTransformBE
    operation QFT(qs : BigEndian) : Unit {
        body (...) {
            ApplyQuantumFourierTransform(BigEndianAsLittleEndian(qs));
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Performs the Quantum Fourier Transform on a quantum register containing an
    /// integer in the little-endian representation.
    ///
    /// # Input
    /// ## qs
    /// Quantum register to which the Quantum Fourier Transform is applied
    ///
    /// # Remarks
    /// The input and output are assumed to be in little endian encoding.
    ///
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.qft"
    operation QFTLE(qs : LittleEndian) : Unit {
        body (...) {
            ApplyQuantumFourierTransform(qs);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Performs a phase shift operation.
    ///
    /// $R=\boldone-(1-e^{i \phi})\ket{1\cdots 1}\bra{1\cdots 1}$.
    ///
    /// # Input
    /// ## phase
    /// The phase $\phi$ applied to state $\ket{1\cdots 1}\bra{1\cdots 1}$.
    /// ## qubits
    /// The register whose state is to be rotated by $R$.
    operation RAll1 (phase : Double, qubits : Qubit[]) : Unit {
        body (...) {
            let nQubits = Length(qubits);
            let flagQubit = qubits[0];
            let systemRegister = qubits[1 .. nQubits - 1];
            Controlled (R1(phase, _))(systemRegister, flagQubit);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Performs a phase shift operation.
    ///
    /// $R=\boldone-(1-e^{i \phi})\ket{0\cdots 0}\bra{0\cdots 0}$.
    ///
    /// # Input
    /// ## phase
    /// The phase $\phi$ applied to state $\ket{0\cdots 0}\bra{0\cdots 0}$.
    /// ## qubits
    /// The register whose state is to be rotated by $R$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RAll1
    operation RAll0 (phase : Double, qubits : Qubit[]) : Unit {
        body (...) {
            ApplyWithCA(ApplyToEachCA(X, _), RAll1(phase, _), qubits);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    /// ## target
    /// Qubit to which the gate should be applied.
    ///
    /// # See Also
    operation HY (target : Qubit) : Unit {
        body (...) {
            H(target);
            S(target);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Applies the Fermionic SWAP.
    ///
    /// # Description
    /// This essentially swaps the qubits while applying a global phase of -1
    /// if both qubits are 1s. Preserves anti-symmetrization of orbitals.
    /// See  for more information.
    ///
    /// This operation is represented by the unitary operator
    /// \begin{align}
    ///     f_{swap} \mathrel{:=}
    ///     \begin{bmatrix}
    ///         1 & 0 & 0 & 0 \\\\
    ///         0 & 0 & 1 & 0 \\\\
    ///         0 & 1 & 0 & 0 \\\\
    ///         0 & 0 & 0 & -1 \\\\
    ///     \end{bmatrix}.
    /// \end{align}
    ///
    /// # Input
    /// ## qubit1
    /// The first qubit to be swapped.
    /// ## qubit2
    /// The second qubit to be swapped.
    ///
    /// # References
    /// - [ *Ryan Babbush, Nathan Wiebe, Jarrod McClean, James McClain,
    ///     Hartmut Neven, Garnet Kin-Lic Chan*,
    ///     arXiv:1706.00023 ](https://arxiv.org/pdf/1706.00023.pdf)
    ///
    /// # See Also
    /// - Microsoft.Quantum.Intrinsic.SWAP
    operation ApplyFermionicSWAP (qubit1 : Qubit, qubit2 : Qubit) : Unit
    is Adj + Ctl {
        SWAP(qubit1, qubit2);
        CZ(qubit1, qubit2);
    }

    /// # Summary
    /// Permutes qubits by using the SWAP operation.
    ///
    /// # Input
    /// ## ordering
    /// Describes the new ordering of the qubits, where the qubit at index i will now be at ordering[i].
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # Example
    /// Given ordering = [2, 1, 0] and register $\ket{\alpha_0} \ket{\alpha_1} \ket{\alpha_2}$, PermuteQubits
    /// changes the register into $\ket{\alpha_2} \ket{\alpha_1} \ket{\alpha_0}$
    ///
    /// ```qsharp
    /// // The following two lines are equivalent
    /// PermuteQubits([2, 1, 0], register);
    /// SWAP(register[0], register[2]);
    /// ```
    operation PermuteQubits(ordering : Int[], register : Qubit[]) : Unit is Adj+Ctl {
        EqualityFactI(Length(ordering), Length(register), "The new ordering has an incorrect number of elements");

        for ((left, right) in SwapOrderToPermuteArray(ordering)) {
            SWAP(register[left], register[right]);
        }
    }

}
