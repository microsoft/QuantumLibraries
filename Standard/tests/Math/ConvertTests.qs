// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;

    @Test("QuantumSimulator")
    function DoubleAsComplexPolarTest() : Unit {
        EqualityFactCP(ComplexPolar(1.0, 0.0), DoubleAsComplexPolar(1.0), "Expected 1.0 == 1.0 exp(0).");
        EqualityFactCP(ComplexPolar(3.0, PI()), DoubleAsComplexPolar(-3.0), "Expected -3.0 == 3.0 exp(𝑖π).");
    }
}
