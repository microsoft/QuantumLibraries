// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Random;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Logical;

    @Test("QuantumSimulator")
    function IntsAreSorted() : Unit {
        Fact(IsSorted(LessThanOrEqualI, [1, 10, 100]), "[1, 10, 100] was marked as unsorted.");
    }

    @Test("QuantumSimulator")
    function IntsAreNotSorted() : Unit {
        Contradiction(IsSorted(LessThanOrEqualI, [100, 10, 3]), "[100, 10, 3] was marked as sorted.");
    }

    @Test("QuantumSimulator")
    function SortedIntsAreSorted() : Unit {
        Fact(IsSorted(LessThanOrEqualI, 
            Sorted(LessThanOrEqualI, [100, 10, 3])),
            "Sorted(<=, [100, 10, 3]) was marked as unsorted."
        );
    }

    @Test("QuantumSimulator")
    operation CheckRandomArraysAreSortedWhenSorted() : Unit {
        let nItems = 100;
        let nTrials = 10;
        let maxItem = 1000;
        for (_ in 0..nTrials - 1) {
            let data = DrawMany((DiscreteUniformDistribution(0, maxItem))::Sample, nItems, ());
            Fact(IsSorted(LessThanOrEqualI, Sorted(LessThanOrEqualI, data)), $"{data} was not sorted after running Sorted.");
        }
    }

}
