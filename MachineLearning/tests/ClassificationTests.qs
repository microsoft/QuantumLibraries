// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.MachineLearning as ML;

    @Test("QuantumSimulator")
    function InferredLabelFact() : Unit {
        EqualityFactI(ML.InferredLabel(0.25, 0.26), 1, "InferredLabel returned wrong class.");
    }

    @Test("QuantumSimulator")
    function InferredLabelsFact() : Unit {
        AllEqualityFactI(
            ML.InferredLabels(0.25, [0.23, 0.26, 0.1]),
            [0, 1, 0],
            "InferredLabels returned at least one wrong class."
        );
    }

}
