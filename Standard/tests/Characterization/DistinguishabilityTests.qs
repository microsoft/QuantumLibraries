// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Diagnostics;

    @Test("QuantumSimulator")
    operation CheckOverlapBetweenPlusAndOne() : Unit {
        let prep1 = ApplyToEachCA(H, _);
        let prep2 = ApplyToEachCA(X, _);

        EqualityWithinToleranceFact(
            EstimateRealOverlapBetweenStates(
                NoOp<Qubit[]>, prep1, prep2, 1, 1000000
            ),
            1.0 / Sqrt(2.0),
            0.02
        );
        EqualityWithinToleranceFact(
            EstimateImagOverlapBetweenStates(
                NoOp<Qubit[]>, prep1, prep2, 1, 1000000
            ),
            0.0,
            0.02
        );
        EqualityWithinToleranceFact(
            0.5,
            EstimateOverlapBetweenStates(
                prep1, prep2, 1, 1000000
            ),
            0.02
        );
    }

    @Test("QuantumSimulator")
    operation CheckOverlapWithCommonPreparation() : Unit {
        let common = ApplyToEachCA(H, _);
        let prep1 = ApplyToEachCA(S, _);
        let prep2 = ApplyToEachCA(Z, _);

        EqualityWithinToleranceFact(
            EstimateRealOverlapBetweenStates(
                common, prep1, prep2, 1, 1000000
            ),
            0.5,
            0.02
        );
        EqualityWithinToleranceFact(
            EstimateImagOverlapBetweenStates(
                common, prep1, prep2, 1, 1000000
            ),
            0.5,
            0.02
        );
    }
}
