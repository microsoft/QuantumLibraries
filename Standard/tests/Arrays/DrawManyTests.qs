// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Random;

    @Test("QuantumSimulator")
    operation CheckDrawMany() : Unit {
        let samples = DrawMany(DrawRandomInt, 20, (0, 29));
        EqualityFactI(20, Length(samples), "Wrong number of samples returned.");

        for sample in samples {
            Fact(0 <= sample and sample < 30, "Sample returned by DrawMany was out of range.");
        }
    }

}
