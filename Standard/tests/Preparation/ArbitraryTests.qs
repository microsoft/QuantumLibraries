// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Preparation.Tests {
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;

    @Test("QuantumSimulator")
    function BlochSphereCoordinatesFact() : Unit {
        let coefficients = (
            ComplexPolar(Sqrt(2.0) / Sqrt(3.0), 0.0),
            ComplexPolar(1.0 / Sqrt(3.0), PI() / 3.0)
        );
        let (prefactor, phi, theta) = BlochSphereCoordinates(coefficients);
        NearEqualityFactCP(prefactor, ComplexPolar(1.0, 0.5235987755982988));
        NearEqualityFactD(phi, 1.0471975511965976);
        NearEqualityFactD(theta, 1.9106332362490186);
    }
}
