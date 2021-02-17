// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Synthesis;

    operation AssertInverseIsCorrect(e : Int, s : Int, x : Int, w : Int) : Unit {
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

    operation AssertOpTimesInverseIsIdentity(e : Int, s : Int, x : Int, w : Int) : Unit {
        Message($"{e} {s} {x} {w}");
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
}
