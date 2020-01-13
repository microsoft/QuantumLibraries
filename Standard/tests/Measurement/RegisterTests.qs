// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Measurement.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Diagnostics;

    @Test("QuantumSimulator")
    operation CheckMeasureIfAllZeros() : Unit {
        using (qs = Qubit[3]) {
            Fact(MeasureIfAllZeros(qs));

            X(qs[1]);
            Fact(not MeasureIfAllZeros(qs));

            ResetAll(qs);
        }
    }

}
