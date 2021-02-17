// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    // Here, we define the group product between different members of the
    // single-qubit Clifford group, using the (𝐸, 𝑆, 𝑋, ω) presentation first
    // used in the newsynth presentation.
    //
    // See the derivation of each commutator and inverse table
    // at /design/notes/sq-clifford-derivations.ipynb.

    internal function XSCommutation(left : SingleQubitClifford, right : SingleQubitClifford) : SingleQubitClifford {
        let xpp = (CanonicalForm1C(left))::X;
        let sp  = (CanonicalForm1C(right))::S;

        if   xpp == 0 and sp == 0 { return SingleQubitClifford(0, 0, 0, 0); }
        elif xpp == 0 and sp == 1 { return SingleQubitClifford(0, 1, 0, 0); }
        elif xpp == 0 and sp == 2 { return SingleQubitClifford(0, 2, 0, 0); }
        elif xpp == 0 and sp == 3 { return SingleQubitClifford(0, 3, 0, 0); }
        elif xpp == 1 and sp == 0 { return SingleQubitClifford(0, 0, 1, 0); }
        elif xpp == 1 and sp == 1 { return SingleQubitClifford(0, 3, 1, 2); }
        elif xpp == 1 and sp == 2 { return SingleQubitClifford(0, 2, 1, 4); }
        elif xpp == 1 and sp == 3 { return SingleQubitClifford(0, 1, 1, 6); }
        else {
            // We guaranteed above that this won't happen, but we can't prove
            // that to the compiler, so we need an else here to catch that
            // case.
            fail "";
        }
    }

    internal function XSECommutation(left : SingleQubitClifford, right : SingleQubitClifford) : SingleQubitClifford {
        let leftCanonical = CanonicalForm1C(left);
        let s = leftCanonical::S;
        let x = leftCanonical::X;
        let ep = (CanonicalForm1C(right))::E;

        if   s == 0 and x == 0 and ep == 0 { return SingleQubitClifford(0, 0, 0, 0); }
        elif s == 0 and x == 0 and ep == 1 { return SingleQubitClifford(1, 0, 0, 0); }
        elif s == 0 and x == 0 and ep == 2 { return SingleQubitClifford(2, 0, 0, 0); }
        elif s == 0 and x == 1 and ep == 0 { return SingleQubitClifford(0, 0, 1, 0); }
        elif s == 0 and x == 1 and ep == 1 { return SingleQubitClifford(1, 2, 0, 0); }
        elif s == 0 and x == 1 and ep == 2 { return SingleQubitClifford(2, 2, 1, 6); }
        elif s == 1 and x == 0 and ep == 0 { return SingleQubitClifford(0, 1, 0, 0); }
        elif s == 1 and x == 0 and ep == 1 { return SingleQubitClifford(2, 3, 0, 6); }
        elif s == 1 and x == 0 and ep == 2 { return SingleQubitClifford(1, 1, 1, 2); }
        elif s == 1 and x == 1 and ep == 0 { return SingleQubitClifford(0, 1, 1, 0); }
        elif s == 1 and x == 1 and ep == 1 { return SingleQubitClifford(2, 1, 0, 6); }
        elif s == 1 and x == 1 and ep == 2 { return SingleQubitClifford(1, 3, 0, 4); }
        elif s == 2 and x == 0 and ep == 0 { return SingleQubitClifford(0, 2, 0, 0); }
        elif s == 2 and x == 0 and ep == 1 { return SingleQubitClifford(1, 2, 1, 6); }
        elif s == 2 and x == 0 and ep == 2 { return SingleQubitClifford(2, 0, 1, 0); }
        elif s == 2 and x == 1 and ep == 0 { return SingleQubitClifford(0, 2, 1, 0); }
        elif s == 2 and x == 1 and ep == 1 { return SingleQubitClifford(1, 0, 1, 2); }
        elif s == 2 and x == 1 and ep == 2 { return SingleQubitClifford(2, 2, 0, 2); }
        elif s == 3 and x == 0 and ep == 0 { return SingleQubitClifford(0, 3, 0, 0); }
        elif s == 3 and x == 0 and ep == 1 { return SingleQubitClifford(2, 1, 1, 4); }
        elif s == 3 and x == 0 and ep == 2 { return SingleQubitClifford(1, 1, 0, 2); }
        elif s == 3 and x == 1 and ep == 0 { return SingleQubitClifford(0, 3, 1, 0); }
        elif s == 3 and x == 1 and ep == 1 { return SingleQubitClifford(2, 3, 1, 0); }
        elif s == 3 and x == 1 and ep == 2 { return SingleQubitClifford(1, 3, 1, 0); }
        else {
            fail "";
        }
    }

    internal function InverseWithoutPhase(op : SingleQubitClifford) : SingleQubitClifford {
        let canonical = CanonicalForm1C(op);
        let e = canonical::E;
        let s = canonical::S;
        let x = canonical::X;

        if   e == 0 and s == 0 and x == 0 { return SingleQubitClifford(0, 0, 0, 0); }
        elif e == 0 and s == 0 and x == 1 { return SingleQubitClifford(0, 0, 1, 0); }
        elif e == 0 and s == 1 and x == 0 { return SingleQubitClifford(0, 3, 0, 0); }
        elif e == 0 and s == 1 and x == 1 { return SingleQubitClifford(0, 1, 1, 6); }
        elif e == 0 and s == 2 and x == 0 { return SingleQubitClifford(0, 2, 0, 0); }
        elif e == 0 and s == 2 and x == 1 { return SingleQubitClifford(0, 2, 1, 4); }
        elif e == 0 and s == 3 and x == 0 { return SingleQubitClifford(0, 1, 0, 0); }
        elif e == 0 and s == 3 and x == 1 { return SingleQubitClifford(0, 3, 1, 2); }
        elif e == 1 and s == 0 and x == 0 { return SingleQubitClifford(2, 0, 0, 0); }
        elif e == 1 and s == 0 and x == 1 { return SingleQubitClifford(2, 2, 1, 6); }
        elif e == 1 and s == 1 and x == 0 { return SingleQubitClifford(1, 1, 0, 2); }
        elif e == 1 and s == 1 and x == 1 { return SingleQubitClifford(1, 3, 0, 2); }
        elif e == 1 and s == 2 and x == 0 { return SingleQubitClifford(2, 0, 1, 0); }
        elif e == 1 and s == 2 and x == 1 { return SingleQubitClifford(2, 2, 0, 6); }
        elif e == 1 and s == 3 and x == 0 { return SingleQubitClifford(1, 1, 1, 2); }
        elif e == 1 and s == 3 and x == 1 { return SingleQubitClifford(1, 3, 1, 2); }
        elif e == 2 and s == 0 and x == 0 { return SingleQubitClifford(1, 0, 0, 0); }
        elif e == 2 and s == 0 and x == 1 { return SingleQubitClifford(1, 2, 0, 0); }
        elif e == 2 and s == 1 and x == 0 { return SingleQubitClifford(2, 1, 1, 4); }
        elif e == 2 and s == 1 and x == 1 { return SingleQubitClifford(2, 1, 0, 4); }
        elif e == 2 and s == 2 and x == 0 { return SingleQubitClifford(1, 2, 1, 6); }
        elif e == 2 and s == 2 and x == 1 { return SingleQubitClifford(1, 0, 1, 6); }
        elif e == 2 and s == 3 and x == 0 { return SingleQubitClifford(2, 3, 0, 6); }
        elif e == 2 and s == 3 and x == 1 { return SingleQubitClifford(2, 3, 1, 2); }
        else {
            fail "";
        }
    }

    /// # Summary
    /// Returns the product of two single-qubit Clifford operators.
    ///
    /// # Input
    /// ## left
    /// The first operator to be multiplied.
    /// ## right
    /// The second operator to be multiplied.
    ///
    /// # Output
    /// The product of `left` and `right`.
    ///
    /// # Example
    /// Suppose that `left` and `right` are both single-qubit Clifford
    /// operators.
    /// ```qsharp
    /// let left = DrawRandomSingleQubitClifford();
    /// let right = DrawRandomSingleQubitClifford();
    /// ```
    /// Then, the following two snippets are equivalent:
    /// ```qsharp
    /// Apply1C(right, q);
    /// Apply1C(left, q);
    /// ```
    /// and:
    /// ```qsharp
    /// Apply1C(Times1C(left, right), q);
    /// ```
    function Times1C(left : SingleQubitClifford, right : SingleQubitClifford) : SingleQubitClifford {
        // Start by finding a new single-qubit Clifford operator
        // 𝑈 = 𝐸^{e''} 𝑆^{s''} 𝑋^{x''} ω^{w''} such that
        // 𝑈 = 𝑆^{s} 𝑋^{x} 𝐸^{e'}, where left! = (e, s, x, w) and where
        // right! = (e', s', x', w'). That is, 𝑈 is the commutation of right's
        // 𝐸 past left's 𝑆 and 𝑋.
        let step1 = XSECommutation(left, right);
        let step2 = XSCommutation(step1, right);

        return CanonicalForm1C(SingleQubitClifford(
            // Combine left's 𝐸 with the contribution we got from commuting
            // right's 𝐸 past left's 𝑆𝑋.
            left::E + step1::E,
            // Next, combine the 𝑆 we got from steps 2 and 3.
            step1::S + step2::S,
            step2::X + right::X,
            step1::Omega + step2::Omega + left::Omega + right::Omega
        ));
    }

    /// # Summary
    /// Returns the inverse of a single-qubit Clifford operators.
    ///
    /// # Input
    /// ## op
    /// The operator to be inverted.
    ///
    /// # Output
    /// The inverse of `op`.
    ///
    /// # Example
    /// Suppose that `op` is a single-qubit Clifford operator.
    /// ```qsharp
    /// let op = DrawRandomSingleQubitClifford();
    /// ```
    /// Then, the following snippet applies the identity (aka no-op) operation:
    /// ```qsharp
    /// Apply1C(op, q);
    /// Apply1C(Inverse1C(op), q);
    /// ```
    function Inverse1C(op : SingleQubitClifford) : SingleQubitClifford {
        let inv = InverseWithoutPhase(op);
        return CanonicalForm1C(inv w/ Omega <- inv::Omega - op::Omega);
    }

    /// # Summary
    /// Returns the action by conjugation of a single-qubit Clifford operator
    /// on a single-qubit Pauli operator.
    ///
    /// # Input
    /// ## op
    /// The single-qubit Clifford operator to conjugate by.
    /// ## pauli
    /// 
    ///
    /// # Example
    /// Given `op` and `pauli`, the following are equivalent up to phase:
    /// ```qsharp
    /// PauliAsSingleQubitClifford(Action1C(op, pauli));
    /// Times1C(op, Times1C(PauliAsSingleQubitClifford(pauli), Inverse1C(op)));
    /// ```
    function Action1C(op : SingleQubitClifford, pauli : Pauli) : Pauli {
        if pauli == PauliI { return PauliI; }

        // Since we're ignoring phase, the X and Omega parts don't do anything.
        // We focus instead on S and E.
        //
        //     𝑆𝑋𝑆⁺ =  𝑌     𝐸𝑋𝐸⁺ =  𝑌
        //     𝑆𝑌𝑆⁺ = −𝑋     𝐸𝑌𝐸⁺ =  𝑍
        //     𝑆𝑍𝑆⁺ =  𝑍     𝐸𝑋𝐸⁺ =  𝑋
        //
        mutable acc = pauli;
        for _ in 0..op::S - 1 {
            if   acc == PauliX { set acc = PauliY; }
            elif acc == PauliY { set acc = PauliX; }
        }
        for _ in 0..op::E - 1 {
            if   acc == PauliX { set acc = PauliY; }
            elif acc == PauliY { set acc = PauliZ; }
            elif acc == PauliZ { set acc = PauliX; }
        }

        return acc;
    }

}
