// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {

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
            RingRepresentative(op::E, 3),
            RingRepresentative(op::S, 4),
            RingRepresentative(op::X, 2),
            RingRepresentative(op::Omega, 8)
        );
    }

    // TODO
    function Identity1C() : SingleQubitClifford {
        return SingleQubitClifford(0, 0, 0, 0);
    }

    // TODO
    function SingleQubitCliffordAsOperation(clifford : SingleQubitClifford) : (Qubit => Unit is Adj + Ctl) {
        return Apply1C(clifford, _);
    }

    // TODO
    function PauliAsSingleQubitClifford(pauli : Pauli) : SingleQubitClifford {
        fail "TODO";
    }

}
