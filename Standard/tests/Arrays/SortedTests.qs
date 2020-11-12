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
    function DoublesAreSorted() : Unit {
        Fact(IsSorted(LessThanOrEqualD, [1.0, 10.1, 100.2]), "[1.0, 10.1, 100.2] was marked as unsorted.");
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
    function SortedDoublesAreSorted() : Unit {
        Fact(IsSorted(LessThanOrEqualD, 
            Sorted(LessThanOrEqualD, [100.0, 10.1, 3.14])),
            "Sorted(<=, [100.0, 10.1, 3.14]) was marked as unsorted."
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

    @Test("QuantumSimulator")
    function LexographicSortIsCorrect() : Unit {
        let arrayComparison = LexographicComparison(LessThanOrEqualI);
        let data = [
            [1, 2, 3],
            [1, 2],
            [0, 2],
            [1, 3]
        ];
        let sorted = Sorted(arrayComparison, data);

        AllEqualityFactI(
            sorted[0], [0, 2], "0th item was not correct."
        );
        AllEqualityFactI(
            sorted[1], [1, 2], "1st item was not correct."
        );
        AllEqualityFactI(
            sorted[2], [1, 2, 3], "2nd item was not correct."
        );
        AllEqualityFactI(
            sorted[3], [1, 3], "3rd item was not correct."
        );
    }

}
