// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Logical;

    @Test("QuantumSimulator")
    operation UniqueInts() : Unit {
        AllEqualityFactI(Unique(EqualI, [0, 0, 1, 2, 2, 3, 4, 5, 5, 8, 42, 42, 39]), [0, 1, 2, 3, 4, 5, 8, 42, 39], "Data is not unique");
        AllEqualityFactI(Unique(EqualI, [0, 1, 1, 0, 0, 1, 1, 0]), [0, 1, 0, 1, 0], "Data is not unique");
        AllEqualityFactI(Unique(EqualI, Sorted(LessThanOrEqualI, [2, 2, 1, 1, 2, 2, 1, 1])), [1, 2], "Sorted data is not unique");
    }

    @Test("QuantumSimulator")
    operation UniqueDoubles() : Unit {
        let unique = Unique(EqualD, [1.1, 1.1, 2.2, 2.2, 2.2, 3.3, 0.5, 42.0]);
        EqualityFactI(Length(unique), 5, "Unexpected length of unique data");
        Fact(unique[0] == 1.1, "Unexpected element in unique data");
        Fact(unique[1] == 2.2, "Unexpected element in unique data");
        Fact(unique[2] == 3.3, "Unexpected element in unique data");
        Fact(unique[3] == 0.5, "Unexpected element in unique data");
        Fact(unique[4] == 42.0, "Unexpected element in unique data");
    }
}
