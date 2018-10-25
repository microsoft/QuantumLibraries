// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Applies the controlled-$X$ (CX) gate to a pair of qubits,
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
    /// ```Q#
    /// Controlled X([control], target);
    /// ```
    ///
    /// Equivalent to:
    /// ```Q#
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
    /// Applies the controlled-$Y$ (CY) gate to a pair of qubits,
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
    /// ```Q#
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
    /// Applies the controlled-$Z$ (CZ) gate to a pair of qubits,
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
    /// ```Q#
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

}