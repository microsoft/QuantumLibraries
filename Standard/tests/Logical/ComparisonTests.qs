// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

    @Test("QuantumSimulator")
    function LexographicComparisonIsCorrect() : Unit {
        let lexographicComparison = LexographicComparison(LessThanOrEqualD);
        Fact(
            lexographicComparison(
                [1.1, 2.2], [1.1, 2.2, 3.3]
            ),
            "Shorter array should have occurred first."
        );
        Fact(
            lexographicComparison(
                [0.7, 2.2], [1.1, 2.2]
            ),
            "Array with smaller first element should have occurred first."
        );
        Fact(
            lexographicComparison(
                [1.1, 2.2], [1.1, 2.2]
            ),
            "Identical arrays should be marked as less than or equal."
        );
        Contradiction(
            lexographicComparison(
                [1.1, 2.7], [1.1, 2.2, 3.3]
            ),
            "Array with larger second element should have occurred second."
        );
    }

}
