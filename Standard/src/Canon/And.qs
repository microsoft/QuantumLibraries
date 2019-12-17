// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Inverts a given target qubit if and only if both control qubits are in the 1 state,
    /// using measurement to perform the adjoint operation.
    ///
    /// # Description
    /// Inverts `target` if and only if both controls are 1, but assumes that
    /// `target` is in state 0.  The operation has T-count 4, T-depth 2 and
    /// requires no helper qubit, and may therefore be preferable to a CCNOT
    /// operation, if `target` is known to be 0.  The adjoint of this operation
    /// is measurement based and requires no T gates.
    ///
    /// # Inputs
    /// ## control1
    /// First control qubit
    /// ## control2
    /// Second control qubit
    /// ## target
    /// Target ancilla qubit; must be in state 0
    ///
    /// # References
    /// - Craig Gidney: "Halving the cost of quantum addition", Quantum 2, page
    ///   74, 2018
    ///   [arXiv:1709.06648](https://arxiv.org/abs/1709.06648)
    ///   doi:10.1103/PhysRevA.85.044302
    operation ApplyAnd(control1 : Qubit, control2 : Qubit, target : Qubit) : Unit {
        body (...) {
            AssertAllZero([target]);
            H(target);
            T(target);
            CNOT(control1, target);
            CNOT(control2, target);
            within {
                CNOT(target, control1);
                CNOT(target, control2);
            }
            apply {
                Adjoint T(control1);
                Adjoint T(control2);
                T(target);
            }
            HY(target);
        }
        adjoint (...) {
            H(target);
            AssertProb([PauliZ], [target], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);
            if (IsResultOne(MResetZ(target))) {
                CZ(control1, control2);
            }
        }
    }

    /// # Summary
    /// Inverts `target` if and only if both controls are 1, but assumes that
    /// `target` is in state 0.  The operation has T-count 4, T-depth 1 and
    /// requires one helper qubit, and may therefore be preferable to a CCNOT
    /// operation, if `target` is known to be 0.  The adjoint of this operation
    /// is measurement based and requires no T gates, and no helper qubit.
    ///
    /// # Inputs
    /// ## control1
    /// First control qubit
    /// ## control2
    /// Second control qubit
    /// ## target
    /// Target ancilla qubit; must be in state 0
    ///
    /// # References
    /// - Craig Gidney: "Halving the cost of quantum addition", Quantum 2, page
    ///   74, 2018
    ///   [arXiv:1709.06648](https://arxiv.org/abs/1709.06648)
    ///   doi:10.1103/PhysRevA.85.044302
    /// - Cody Jones: "Novel constructions for the fault-tolerant Toffoli gate",
    ///   Phys. Rev. A 87, 022328, 2013
    ///   [arXiv:1212.5069](https://arxiv.org/abs/1212.5069)
    ///   doi:10.1103/PhysRevA.87.022328
    operation ApplyLowDepthAnd(control1 : Qubit, control2 : Qubit, target : Qubit) : Unit {
        body (...) {
            using (helper = Qubit()) {
                AssertAllZero([target]);
                H(target);
                within {
                    CNOT(control1, helper);
                    CNOT(control2, helper);
                    CNOT(target, control1);
                    CNOT(target, control2);
                }
                apply {
                    Adjoint T(control1);
                    Adjoint T(control2);
                    T(target);
                    T(helper);
                }
                HY(target);
            }
        }
        adjoint (...) {
            Adjoint ApplyAnd(control1, control2, target);
        }
    }
}
