// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.MachineLearning as ML;

    @Test("QuantumSimulator")
    function MisclassificationsFact() : Unit {
        let misclassifications = ML.Misclassifications([0, 1, 0, 0], [0, 1, 1, 0]);
        AllEqualityFactI(misclassifications, [2], "Wrong output from Misclassifications.");
    }

    @Test("QuantumSimulator")
    function NMisclassificationsFact() : Unit {
        let nMisclassifications = ML.NMisclassifications([0, 1, 0, 0, 1], [0, 1, 1, 0, 0]);
        EqualityFactI(nMisclassifications, 2, "Wrong output from NMisclassifications.");
    }
}
