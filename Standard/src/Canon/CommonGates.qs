// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Arithmetic;

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
            ApplyQuantumFourierTransformBE(qs);
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
            ApplyQuantumFourierTransformLE(qs);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

}
