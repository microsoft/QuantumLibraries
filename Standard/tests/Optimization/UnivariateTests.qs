// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Optimization;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;

    function ParabolaCase(minima : Double, x : Double) : Double {
        return PowD((x - minima), 2.0);
    }

    @Test("QuantumSimulator")
    function MinimizedParabolaIsCorrect() : Unit {
        let optimum = LocalUnivariateMinimum(ParabolaCase(3.14, _), (-7.0, +12.0), 1e-10);
        NearEqualityFactD(optimum::Coordinate, 3.14);
        NearEqualityFactD(optimum::Value, 0.0);
    }

}
