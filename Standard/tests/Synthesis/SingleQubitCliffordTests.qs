// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Synthesis;

    internal operation AssertInverseIsCorrect(e : Int, s : Int, x : Int, w : Int) : Unit {
        let c = SingleQubitClifford(e, s, x, w);
        let op = Apply1C(c, _);
        let inv = Inverse1C(c);
        let invOp = Apply1C(inv, _);

        AssertOperationsEqualReferenced(1,
            NoOp<Qubit[]>,
            ApplyToFirstQubitCA(BoundCA([op, invOp]), _)
        );
    }


    /// # Summary
    /// Tests that applying a single qubit Clifford and then applying its inverse
    /// correctly results in no evolution of the target qubit.
    ///
    /// # Remarks
    /// This does not test that the group product on the single-qubit Clifford
    /// group is correct.
    @Test("QuantumSimulator")
    operation AssertAllInversesAreCorrect() : Unit {
        for e in 0..2 {
            for s in 0..3 {
                for x in 0..1 {
                    for w in 0..7 {
                        AssertInverseIsCorrect(e, s, x, w);
                    }
                }
            }
        }
    }

    
    @Test("QuantumSimulator")
    function TimesInverseIsIdentity() : Unit {
        for e in 0..2 {
            for s in 0..3 {
                for x in 0..1 {
                    for w in 0..7 {
                        let c = SingleQubitClifford(e, s, x, w);
                        let invC = Inverse1C(c);
                        IdentityFact1C(Times1C(c, invC), $"c × Adj c was wrong for c = {c}");
                        IdentityFact1C(Times1C(invC, c), $"Adj c × c was wrong for c = {c}");
                    }
                }
            }
        }
    }

    internal operation AssertOpTimesInverseIsIdentity(e : Int, s : Int, x : Int, w : Int) : Unit {
        let c = SingleQubitClifford(e, s, x, w);
        let op = Apply1C(Times1C(c, Inverse1C(c)), _);

        AssertOperationsEqualReferenced(1,
            NoOp<Qubit[]>,
            ApplyToFirstQubitCA(op, _)
        );
    }

    @Test("QuantumSimulator")
    operation AssertAllOpsTimesInverseAreIdentity() : Unit {
        for e in 0..2 {
            for s in 0..3 {
                for x in 0..1 {
                    for w in 0..7 {
                        AssertOpTimesInverseIsIdentity(e, s, x, w);
                    }
                }
            }
        }
    }

    internal operation AssertProductIsCorrect(e : Int, s : Int, x : Int, ep : Int, sp : Int, xp : Int)
    : Unit {
        let left = SingleQubitClifford(e, s, x, 0);
        let right = SingleQubitClifford(ep, sp, xp, 0);
        let product = Times1C(left, right);
        let productAsOp = Apply1C(Times1C(left, right), _);
        let bound = BoundCA(Mapped(SingleQubitCliffordAsOperation, [right, left]));

        AssertOperationsEqualReferenced(1,
            ApplyToFirstQubitCA(bound, _),
            ApplyToFirstQubitCA(productAsOp, _)
        );
    }

    @Test("QuantumSimulator")
    operation AssertAllProductsWithoutPhasesAreCorrect() : Unit {
        for e in 0..2 {
            for s in 0..3 {
                for x in 0..1 {
                    for ep in 0..3 {
                        for sp in 0..3 {
                            for xp in 0..1 {
                                AssertProductIsCorrect(e, s, x, ep, sp, xp);
                            }
                        }
                    }
                }
            }
        }
    }
}
