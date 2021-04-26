// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Math;

    /// # Summary
    /// An element of the single-qubit Clifford group.
    ///
    /// # Description
    /// Each element is represented as as $E^{i} S^{j} X^{k} \omega^{\ell}$, where $\omega^8 = 1$ and where $E = H S^3 \omega^3$.
    ///
    /// # Named Items
    /// ## E
    /// The power to which $E$ should be raised to generate this element.
    /// ## S
    /// The power to which $S$ should be raised to generate this element.
    /// ## X
    /// The power to which $X$ should be raised to generate this element.
    /// ## Omega
    /// The power to which $\omega$ should be raised to generate this element.
    ///
    /// # Example
    /// ```Q#
    /// let identity = SingleQubitClifford((0, 0, 0, 0));
    /// let xClifford = SingleQubitClifford((0, 0, 1, 0));
    /// ```
    ///
    /// # References
    /// - https://hackage.haskell.org/package/newsynth-0.4.0.0/docs/Quantum-Synthesis-Clifford.html
    newtype SingleQubitClifford = (
        E: Int,
        S: Int,
        X: Int,
        Omega: Int
    );

    /// # Summary
    /// Returns a representation of a given single-qubit Clifford operator
    /// in its canonical form.
    ///
    /// # Description
    /// The canonical form is defined such that if two values of type
    /// `SingleQubitClifford` have the same canonical form if and only if they
    /// represent the same single-qubit Clifford operator.
    internal function CanonicalForm1C(op : SingleQubitClifford) : SingleQubitClifford {
        return SingleQubitClifford(
            ModulusI(op::E, 3),
            ModulusI(op::S, 4),
            ModulusI(op::X, 2),
            ModulusI(op::Omega, 8)
        );
    }

    /// # Summary
    /// Returns a representation of the identity as a single-qubit Clifford
    /// operator.
    function Identity1C() : SingleQubitClifford {
        return SingleQubitClifford(0, 0, 0, 0);
    }

    /// # Summary
    /// Returns a representation of a single-qubit Clifford operator as an
    /// operation acting on a single qubit.
    ///
    /// # Input
    /// The operator to be represented as an operation.
    ///
    /// # Output
    /// An operation that applies the given Clifford operator to a single
    /// qubit.
    ///
    /// # Example
    /// Suppose that `op` is a single-qubit Clifford operator, and that
    /// `q` is a single qubit:
    ///
    /// ```qsharp
    /// let op = DrawRandomSingleQubitClifford();
    /// use q = Qubit();
    /// ```
    ///
    /// Then, the following two lines are equivalent:
    /// ```qsharp
    /// Apply1C(op, q);
    /// SingleQubitCliffordAsOperation(op)(q);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.Apply1C
    function SingleQubitCliffordAsOperation(clifford : SingleQubitClifford) : (Qubit => Unit is Adj + Ctl) {
        return Apply1C(clifford, _);
    }

    /// # Summary
    /// Returns a representation of a single-qubit Pauli operator as
    /// a single-qubit Clifford operator that acts by conjugation.
    ///
    /// # Description
    /// Given a Pauli operator $P$, this function returns a Clifford operator
    /// that represents the function $Q \mapsto PQP^{\dagger}$.
    ///
    /// # Input
    /// ## pauli
    /// The Pauli operator to be represented as a Clifford operator.
    ///
    /// # Ouput
    /// A single-qubit Clifford operator representing the action of `pauli` by
    /// conjugation.
    ///
    /// # Remarks
    /// For a value `pauli` of type `Pauli`, `ApplyP(pauli, _)` and
    /// `Apply1C(PauliAsSingleQubitClifford(pauli), _)` are the same operation.
    function PauliAsSingleQubitClifford(pauli : Pauli) : SingleQubitClifford {
        if   pauli == PauliI { return Identity1C(); }
        elif pauli == PauliX { return Identity1C() w/ X <- 1; }
        elif pauli == PauliY { return Identity1C() w/ X <- 1 w/ Omega <- 6 w/ S <- 2; }
        elif pauli == PauliZ { return Identity1C() w/ S <- 2; }
        else {
            fail "";
        }
    }

}
