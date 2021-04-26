// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Optimization;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;

    function ParabolaCase(minima : Double, x : Double) : Double {
        return PowD((x - minima), 2.0);
    }
    
    function AbsDistance(x : Double, y : Double) : Double {
        return AbsD(x - y);
    }

    @Test("QuantumSimulator")
    function MinimizedParabolaIsCorrect() : Unit {
        let optimum = LocalUnivariateMinimum(ParabolaCase(3.14, _), (-7.0, +12.0), 1e-10);
        NearEqualityFactD(optimum::Coordinate, 3.14);
        NearEqualityFactD(optimum::Value, 0.0);
        EqualityFactI(optimum::NQueries, 56, "LocalUnivariateMinimum made an unexpected amount of queries.");
    }

    @Test("QuantumSimulator")
    function MinimumConvergesToEdge() : Unit {
        let optimum = LocalUnivariateMinimum(AbsDistance(-1.0, _), (-1.0, 1.0), 1e-2);
        EqualityWithinToleranceFact(optimum::Coordinate, -1.0, 1e-2);
        EqualityFactI(optimum::NQueries, 14, "LocalUnivariateMinimum made an unexpected amount of queries.");
    }
}
