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

    // TODO
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

    // TODO
    function Inverse1C(op : SingleQubitClifford) : SingleQubitClifford {
        let inv = InverseWithoutPhase(op);
        return CanonicalForm1C(inv w/ Omega <- inv::Omega - op::Omega);
    }

    // TODO
    function Action1C(op : SingleQubitClifford, pauli : Pauli) : Pauli {
        fail "TODO";
    }

}
